"""Проверка инструмента через настоящий stdio MCP-сервер (без агента).

Запускает файл запуска агента `bin/family-profile-mcp` как отдельный процесс,
подключается по протоколу MCP, перечисляет инструменты и вызывает `family_profile`.

Запуск (интерпретатором окружения агента, где стоит пакет mcp):
    /Users/elizavetamenshchikova/YSDA/AI/AI4Science/.venv/bin/python \
        hw3/tool/call_family_profile_stdio.py Kr2_S81
"""
from __future__ import annotations

import asyncio
import json
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

LAUNCHER = Path("/Users/elizavetamenshchikova/YSDA/AI/AI4Science/bin/family-profile-mcp")


async def main(sample_id: str) -> int:
    params = StdioServerParameters(command=str(LAUNCHER), args=[])
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            listed = await session.list_tools()
            names = [t.name for t in listed.tools]
            print("инструменты сервера:", names)
            if "family_profile" not in names:
                print("FAIL: family_profile не опубликован")
                return 1
            result = await session.call_tool("family_profile", {"sample_id": sample_id})
            print("isError:", result.isError)
            for block in result.content:
                text = getattr(block, "text", None)
                if text is not None:
                    print(text)
            return 0 if not result.isError else 1


if __name__ == "__main__":
    sid = sys.argv[1] if len(sys.argv) > 1 else "Kr2_S81"
    sys.exit(asyncio.run(main(sid)))
