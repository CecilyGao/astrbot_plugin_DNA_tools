import re
from typing import List, Optional, Tuple
from astrbot.api.event import filter, AstrMessageEvent
from astrbot.api.star import Context, Star, register
from astrbot.api import logger, AstrBotConfig
# ======================== 常量 ========================
ALLOWED_BASES = set("ATCGUN")
# ======================== 主插件类 ========================
@register(
    "astrbot_plugin_DNA_tools",
    "CecilyGao",
    "参考擎科生物官网云工具，实现DNA格式化、反向互补、翻译、gc含量、引物分析、酶切位点分析功能",
    "0.1.0",
    "https://github.com/CecilyGao/astrbot_plugin_DNA_tools"
)
class DNAPlugin(Star):
    def __init__(self, context: Context, config: AstrBotConfig):
        super().__init__(context)
        self.config = config
        # 群白名单配置（留空则所有群可用）
        self.group_whitelist = config.get('group_whitelist', [])
        # 格式化每行宽度（可配置）
        self.format_width = config.get('format_width', 60)

    # ======================== 权限辅助 ========================
    def _check_group_permission(self, event: AstrMessageEvent) -> bool:
        if not self.group_whitelist:
            return True
        group_id = str(event.get_group_id())
        return group_id in self.group_whitelist

    # ======================== 主命令入口（手动解析） ========================
    @filter.command("DNA", alias={'dna'})
    async def dna_command(self, event: AstrMessageEvent):
        if not self._check_group_permission(event):
            yield event.plain_result("❌ 当前群未在DNA工具白名单中。")
            return

        # 获取完整消息并去除命令前缀（例如 "DNA " 或 "dna "）
        full_text = event.message_str.strip()
        # 移除命令本身（不区分大小写），只保留子命令和参数
        # 由于命令可能带有别名，我们通过正则或简单分割来处理
        # 简单方法：按空格分割，第一个是命令，后面的是子命令和参数
        parts = full_text.split()
        if len(parts) < 2:
            # 只有命令，没有子命令 -> 显示帮助
            yield event.plain_result(self._get_help())
            return

        subcmd = parts[1].lower()
        args = parts[2:] if len(parts) > 2 else []

        # 帮助命令
        if subcmd in ['help', '帮助']:
            yield event.plain_result(self._get_help())
            return

        # 需要序列的子命令列表
        seq_required_cmds = ['format', '格式化', 'reverse', '反向互补', 
                             'translate', '翻译', 'gc', 'gc含量',
                             'primer', '引物分析', 'restriction', '酶切']
        if subcmd in seq_required_cmds:
            seq = self._extract_sequence(args)
            if seq is None:
                yield event.plain_result("⚠️ 请提供有效的DNA序列（只允许 A、T、C、G、U、N），例如：`DNA format ATCG`")
                return
        else:
            # 未知子命令
            yield event.plain_result(f"❌ 未知子命令：{subcmd}\n" + self._get_help())
            return

        # 子命令分发
        if subcmd in ['format', '格式化']:
            result = await self._cmd_format(seq)
            yield event.plain_result(result)
        elif subcmd in ['reverse', '反向互补']:
            result = await self._cmd_reverse(seq)
            yield event.plain_result(result)
        elif subcmd in ['translate', '翻译']:
            result = await self._cmd_translate(seq)
            yield event.plain_result(result)
        elif subcmd in ['gc', 'gc含量']:
            result = await self._cmd_gc(seq)
            yield event.plain_result(result)
        elif subcmd in ['primer', '引物分析']:
            result = await self._cmd_primer(seq)
            yield event.plain_result(result)
        elif subcmd in ['restriction', '酶切']:
            result = await self._cmd_restriction(seq)
            yield event.plain_result(result)
        else:
            yield event.plain_result(f"❌ 未知子命令：{subcmd}\n" + self._get_help())

    # ======================== 参数提取辅助 ========================
    @staticmethod
    def _extract_sequence(args: List[str]) -> Optional[str]:
        """从命令参数列表中提取并清洗DNA序列"""
        if not args:
            return None
        raw = ''.join(args).replace(' ', '').replace('\n', '').replace('\r', '')
        seq = raw.upper()
        if not seq:
            return None
        # 校验序列（只允许 A、T、C、G、U、N）
        if not set(seq).issubset(ALLOWED_BASES):
            return None
        return seq

    # ======================== 子命令实现（保持不变） ========================
    async def _cmd_format(self, seq: str) -> str:
        width = self.format_width
        formatted = '\n'.join([seq[i:i+width] for i in range(0, len(seq), width)])
        return f"📄 格式化结果（每行{width}个碱基）：\n{formatted}"

    async def _cmd_reverse(self, seq: str) -> str:
        complement = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'U': 'A', 'N': 'N'
        }
        rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
        return f"🧬 反向互补序列：\n{rev_comp}"

    async def _cmd_translate(self, seq: str) -> str:
        codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        }
        protein = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3].upper()
            protein.append(codon_table.get(codon, 'X'))
        return f"🧪 翻译结果（标准密码子表）：\n{''.join(protein)}"

    async def _cmd_gc(self, seq: str) -> str:
        gc = self._gc_content(seq)
        repeats = self._find_simple_repeats(seq)
        repeat_str = ', '.join(repeats[:5]) if repeats else "无"
        return f"📊 GC含量：{gc:.2f}%\n🔁 重复序列（长度2-6，仅显示前5个）：{repeat_str}"

    async def _cmd_primer(self, seq: str) -> str:
        gc = self._gc_content(seq)
        tm = self._calculate_tm(seq)
        return (f"🧬 引物分析：\n"
                f"序列：{seq}\n"
                f"GC含量：{gc:.2f}%\n"
                f"Tm值（简单公式）：{tm:.2f} °C")

    async def _cmd_restriction(self, seq: str) -> str:
        sites = self._find_restriction_sites(seq)
        if sites:
            lines = ["🔬 发现的酶切位点："]
            for enzyme, site, pos in sites:
                lines.append(f"  {enzyme} 识别序列 {site} 位置 {pos}")
            return "\n".join(lines)
        else:
            return "❌ 未发现已知酶切位点。"

    # ======================== 核心功能函数 ========================
    def _gc_content(self, seq: str) -> float:
        seq = seq.upper()
        if not seq:
            return 0.0
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq) * 100

    def _find_simple_repeats(self, seq: str, max_len: int = 6) -> List[str]:
        repeats = set()
        seq_upper = seq.upper()
        for length in range(2, max_len + 1):
            for i in range(len(seq_upper) - length + 1):
                sub = seq_upper[i:i+length]
                if seq_upper.count(sub) > 1:
                    repeats.add(sub)
        return sorted(list(repeats), key=lambda x: len(x), reverse=True)

    def _calculate_tm(self, seq: str) -> float:
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        at = len(seq) - gc
        return 4 * gc + 2 * at

    def _find_restriction_sites(self, seq: str) -> List[Tuple[str, str, int]]:
        enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'NotI': 'GCGGCCGC',
            'XbaI': 'TCTAGA',
            'SalI': 'GTCGAC',
            'PstI': 'CTGCAG',
            'SmaI': 'CCCGGG',
            'KpnI': 'GGTACC',
            'SacI': 'GAGCTC',
        }
        found = []
        seq_upper = seq.upper()
        for enzyme, site in enzymes.items():
            pos = seq_upper.find(site)
            if pos != -1:
                found.append((enzyme, site, pos))
        return found

    # ======================== 帮助信息 ========================
    def _get_help(self) -> str:
        return (
            "🧬 DNA工具使用帮助\n"
            "命令格式：DNA <子命令> [序列]\n\n"
            "可用子命令（支持英文或中文别名）：\n"
            "  format / 格式化   <序列>  – 每行60个碱基格式化输出\n"
            "  reverse / 反向互补 <序列>  – 计算反向互补序列\n"
            "  translate / 翻译  <序列>  – 使用标准密码子表翻译为氨基酸\n"
            "  gc / gc含量       <序列>  – 计算GC含量并查找简单重复序列\n"
            "  primer / 引物分析 <序列>  – 计算GC%和Tm值（短引物公式）\n"
            "  restriction / 酶切 <序列> – 查找常见酶切位点（EcoRI, BamHI等）\n"
            "  help / 帮助              – 显示此帮助\n\n"
            "示例：\n"
            "  DNA format ATCGATCGATCG\n"
            "  DNA 反向互补 ATGC\n"
            "  DNA translate ATGCCGTAA\n"
            "  DNA gc ATGCGC\n"
            "  DNA primer ATGCGT\n"
            "  DNA 酶切 GAATTC\n\n"
            "注意：序列中的空格和换行会被自动忽略，只允许字母 ATCGUN。"
        )