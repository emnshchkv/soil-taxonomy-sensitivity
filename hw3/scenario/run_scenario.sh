#!/usr/bin/env bash
# Сценарий «генератор — критик» без участия человека.
#
#   bash scenario/run_scenario.sh                         # обычный запуск
#   tmux new -d -s hw3 'bash scenario/run_scenario.sh'    # запуск, который переживает обрыв связи
#   tmux attach -t hw3                                    # вернуться к сеансу
#
# Роли по очереди пишут на общую шину bus/: генератор — draft_N.md, критик — review_N.md.
# Сценарий останавливается, когда критик принял черновик или исчерпаны раунды.
# Состояние хранится в scenario/state: прерванный сценарий при повторном запуске продолжается.
set -euo pipefail
cd "$(dirname "$0")/.."

MODEL="${MODEL:-deepseek}"
MAX_ROUNDS="${MAX_ROUNDS:-3}"
BUS=bus
LOG=scenario/journal.log
STATE=scenario/state
mkdir -p "$BUS"

log() { printf '%s %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"; }

save_state() {
  printf '%s %s\n' "$1" "$2" > "$STATE.tmp"
  mv "$STATE.tmp" "$STATE"
}

# запись на шину регистрируется с контрольной суммой: по ней видно, кто и когда что написал
register() {
  local role="$1" file="$2" sum
  sum="$(shasum -a 256 "$file" | cut -c1-16)"
  printf '%s | %s | %s | sha256 %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$role" "$file" "$sum" >> "$BUS/INDEX.md"
}

role() {   # role <имя> <раунд> <ожидаемый файл>
  local name="$1" round="$2" output="$3" prev=$(( $2 - 1 )) rc=0 saved
  local prompt; prompt="$(sed -e "s/{round}/$round/g" -e "s/{prev}/$prev/g" "scenario/roles/$name.md")"
  if [ -e "$output" ]; then
    saved="$(mktemp -d "scenario/retry_${name}_${round}.XXXXXX")"
    mv "$output" "$saved/"
    log "раунд $round: предыдущий результат $name сохранён в $saved"
  fi
  log "раунд $round: $name начал"
  ai4science run -m "$MODEL" "$prompt" >> "scenario/run_${name}_${round}.log" 2>&1 || rc=$?
  log "раунд $round: $name закончил, код $rc"
  if [ "$rc" -ne 0 ]; then
    log "остановка: при повторном запуске будет повторён шаг $name раунда $round"
    return "$rc"
  fi
  if [ ! -f "$output" ] || [ ! -s "$output" ]; then
    log "ошибка: шаг $name раунда $round не создал непустой файл $output"
    return 1
  fi
}

round=1; step=generator
if [ -f "$STATE" ]; then
  read -r round step < "$STATE"
  if [ "$round" = done ]; then
    log "сценарий уже завершён: $step"
    [[ "$step" == принято* ]] && exit 0
    exit 1
  fi
  log "продолжение после остановки: раунд $round, шаг $step"
else
  log "начало сценария: модель $MODEL, раундов не больше $MAX_ROUNDS"
  save_state "$round" "$step"
fi

while [ "$round" -le "$MAX_ROUNDS" ]; do
  if [ "$step" = generator ]; then
    role generator "$round" "$BUS/draft_$round.md"
    register генератор "$BUS/draft_$round.md"
    step=critic; save_state "$round" "$step"
  fi
  review="$BUS/review_$round.md"
  role critic "$round" "$review"
  verdict="$(head -n 1 "$review")"
  case "$verdict" in
    'ВЕРДИКТ: принято'|'ВЕРДИКТ: отклонено') register критик "$review" ;;
    *) log "ошибка: неверная первая строка $review; при запуске повторится шаг critic"; exit 1 ;;
  esac
  if [ "$verdict" = 'ВЕРДИКТ: принято' ]; then
    log "итог: черновик принят в раунде $round"
    save_state done "принято в раунде $round"
    exit 0
  fi
  log "раунд $round: критик отклонил черновик"
  round=$(( round + 1 )); step=generator; save_state "$round" "$step"
done
log "итог: раунды исчерпаны, черновик не принят"
save_state done "не принят за $MAX_ROUNDS раундов"
exit 1
