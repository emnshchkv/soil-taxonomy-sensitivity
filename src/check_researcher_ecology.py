"""Read-only check of the user-supplied ecology transcription and sample mapping."""
from qiime_runtime import selected_python
selected_python()

import csv
import hashlib
import json
from pathlib import Path
import yaml


def main():
    root = Path(__file__).resolve().parents[1]
    record = yaml.safe_load((root / 'data/meta/researcher_ecology_20261118.yaml').read_text())
    errors = []

    def check(ok, message):
        if not ok:
            errors.append(message)

    expected_map = {f'Kr{i}_S{79+i}': f'КР{i}' for i in range(1, 6)}
    expected_map.update({f'LUB{i}_S{90+i}': 'ЛЮБ1' if i <= 4 else 'ЛЮБ2' for i in range(2, 7)})
    mappings = record['sample_location_map']
    actual = {row['sample_id']: row['location_code'] for row in mappings}
    check(actual == expected_map and len(mappings) == len(expected_map), 'Привязка sample-id не совпадает с ответом пользователя')
    with (root / 'data/meta/metadata.tsv').open() as stream:
        metadata = list(csv.DictReader(stream, delimiter='\t'))
    check(set(actual) == {row['sample-id'] for row in metadata if row['sample_type'] == 'soil'}, 'Не совпадает набор почвенных образцов')
    check(record['sampling_depth_cm'] == 8, 'Глубина должна быть 8 см')
    check(record['negative_control'] == {'sample_id': 'K_S63', 'sampling_depth_cm': None, 'location_code': None}, 'Контролю нельзя приписывать почвенный адрес или глубину')
    expected_addresses = {'КР1': 'Крылатские холмы 45, корпус 2', 'КР2': 'Осенний бульвар 16к1',
        'КР3': 'Осенний бульвар 7к4', 'КР4': 'Рублёвское ш., 48/1', 'КР5': 'Осенний б-р, 5к3',
        'КР6': 'Крылатская ул., д.1 соор.1', 'ЛЮБ1': 'Судакова 29', 'ЛЮБ2': 'Совхозная 41, детская площадка'}
    addresses = {row['code']: row['address'] for row in record['locations']}
    check(addresses == expected_addresses and len(record['locations']) == 8, 'Адреса не совпадают с сообщением')
    check(record['location_without_sequencing_sample'] == ['КР6'], 'КР6 должна оставаться точкой без файла ридов')
    expected_chem = {'КР1': ['0,110', '0,430', '0,067', '5,370'], 'ЛЮБ1': ['0,050', '1,870', '0,086', '2,460']}
    chemical_records = record['chemistry']
    check(len(chemical_records) == 2 and {c['location_code'] for c in chemical_records} == set(expected_chem), 'Химия должна относиться только к КР1 и ЛЮБ1')
    for site in chemical_records:
        code = site['location_code']
        expected = [{'label_on_image': label, 'unit': unit, 'value': value.replace(',', '.'), 'source_display': value}
                    for label, unit, value in zip(['NH4+', 'NO3-', 'NO2-', 'Organic carbon'], ['г/л', 'г/л', 'г/л', '%'], expected_chem.get(code, []))]
        check(site['measurements'] == expected, f'{code}: транскрипция значений или единиц отличается от изображения')
        check(site['extract_type'] is None and site['measurement_method'] is None, f'{code}: неизвестный метод/вытяжка не должны заполняться догадкой')
        check('одна проба' in site['sampling_basis'], f'{code}: не указан индивидуальный характер химической пробы')
    image = root / record['image_file']
    check(image.is_file() and hashlib.sha256(image.read_bytes()).hexdigest() == record['image_sha256'], 'Исходное изображение отсутствует или изменено')
    for sample, code in expected_map.items():
        note = (root / 'notes' / f'Образец {sample}.md').read_text()
        check('Глубина отбора почвы — 8 см.' in note and f'Код точки — {code};' in note and addresses[code] in note,
              f'{sample}: заметка не согласована с дополнительными метаданными')
    result = {'status': 'FAIL' if errors else 'PASS', 'errors': errors, 'mapped_soil_samples': len(actual),
              'locations': len(addresses), 'chemistry_locations': list(expected_chem), 'chemistry_values': 8,
              'scope': 'transcription and user-confirmed mapping; not chemical-method or ecological validation'}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
