import re
from astrbot.api.event import filter, AstrMessageEvent
from astrbot.api.star import Context, Star, register
from astrbot.api import logger
from astrbot.api import AstrBotConfig

@register(
    "astrbot_plugin_DNA_tools",
    "CecilyGao",
    "参考擎科生物官网云工具，实现DNA格式化、反向互补、翻译、gc含量、引物分析、酶切位点分析功能",
    "0.0.1",
    "https://github.com/CecilyGao/astrbot_plugin_DNA_tools"
)
class DNAPlugin(Star):
    def __init__(self, context: Context, config: AstrBotConfig):
        super().__init__(context)
        self.context = context
        self.config = config
        # 群白名单配置（留空则所有群可用）
        self.group_whitelist = config.get('group_whitelist', [])
        # 格式化每行宽度（可配置）
        self.format_width = config.get('format_width', 60)

    def _check_group_permission(self, event: AstrMessageEvent) -> bool:
        if not self.group_whitelist:
            return True
        group_id = str(event.get_group_id())
        return group_id in self.group_whitelist

    # ---------- 主命令 ----------
    @filter.command("DNA", alias={'dna'})
    async def dna_command(self, event: AstrMessageEvent, subcommand: str = None, *args):
        if not self._check_group_permission(event):
            yield event.plain_result("❌ 当前群未在DNA工具白名单中。")
            return

        # 无子命令或 help 则显示帮助
        if not subcommand or subcommand.lower() in ['help', '帮助']:
            yield event.plain_result(self.get_help())
            return

        # 提取序列参数（去除空格和换行，转为大写）
        seq = None
        if args:
            # 将所有参数拼接成一个字符串，并移除空白字符
            seq = ''.join(args).replace(' ', '').replace('\n', '').replace('\r', '').upper()
        if not seq:
            yield event.plain_result("⚠️ 请提供DNA序列作为参数，例如：`DNA format ATCGATCG`\n" + self.get_help())
            return

        # 校验序列（只允许 ATCGUN 等）
        if not re.match(r'^[ATCGUN]+$', seq):
            yield event.plain_result("❌ 序列包含非法字符，只允许 A、T、C、G、U、N（大小写均可）。")
            return

        cmd = subcommand.lower()

        # ---------- 子命令分发 ----------
        if cmd in ['format', '格式化']:
            result = self.format_sequence(seq)
            yield event.plain_result(f"📄 格式化结果（每行{self.format_width}个碱基）：\n{result}")

        elif cmd in ['reverse', '反向互补']:
            result = self.reverse_complement(seq)
            yield event.plain_result(f"🧬 反向互补序列：\n{result}")

        elif cmd in ['translate', '翻译']:
            protein = self.translate_dna(seq)
            yield event.plain_result(f"🧪 翻译结果（标准密码子表）：\n{protein}")

        elif cmd in ['gc', 'gc含量']:
            gc = self.gc_content(seq)
            repeats = self.find_simple_repeats(seq)
            repeat_str = ', '.join(repeats[:5]) if repeats else "无"
            yield event.plain_result(f"📊 GC含量：{gc:.2f}%\n🔁 重复序列（长度2-6，仅显示前5个）：{repeat_str}")

        elif cmd in ['primer', '引物分析']:
            tm = self.calculate_tm(seq)
            gc = self.gc_content(seq)
            yield event.plain_result(f"🧬 引物分析：\n序列：{seq}\nGC含量：{gc:.2f}%\nTm值（简单公式）：{tm:.2f} °C")

        elif cmd in ['restriction', '酶切']:
            sites = self.find_restriction_sites(seq)
            if sites:
                lines = ["🔬 发现的酶切位点："]
                for enzyme, site, pos in sites:
                    lines.append(f"  {enzyme} 识别序列 {site} 位置 {pos}")
                yield event.plain_result("\n".join(lines))
            else:
                yield event.plain_result("❌ 未发现已知酶切位点。")

        else:
            yield event.plain_result(f"❌ 未知子命令：{subcommand}\n" + self.get_help())

    # ---------- 帮助信息 ----------
    def get_help(self) -> str:
        return """
🧬 DNA工具使用帮助
母命令：DNA（或 dna）
子命令（支持英文或中文别名）：
  format / 格式化   <序列>   – 每行60个碱基格式化输出
  reverse / 反向互补 <序列>  – 计算反向互补序列
  translate / 翻译  <序列>   – 使用标准密码子表翻译为氨基酸
  gc / gc含量       <序列>   – 计算GC含量并查找简单重复序列
  primer / 引物分析 <序列>   – 计算GC%和Tm值（短引物公式）
  restriction / 酶切 <序列>  – 查找常见酶切位点（EcoRI, BamHI等）
  help / 帮助                – 显示此帮助

示例：
  DNA format ATCGATCGATCG
  DNA 反向互补 ATGC
  DNA translate ATGCCGTAA
  DNA gc ATGCGC
  DNA primer ATGCGT
  DNA 酶切 GAATTC
注意：序列中的空格和换行会被自动忽略，只允许字母 ATCGUN。
"""

    # ---------- 核心功能函数 ----------
    def format_sequence(self, seq: str, width: int = None) -> str:
        if width is None:
            width = self.format_width
        return '\n'.join([seq[i:i+width] for i in range(0, len(seq), width)])

    def reverse_complement(self, seq: str) -> str:
        complement = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'U': 'A', 'N': 'N'
        }
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def translate_dna(self, seq: str) -> str:
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
        # 从第一个碱基开始，每次取3个
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3].upper()
            protein.append(codon_table.get(codon, 'X'))
        return ''.join(protein)

    def gc_content(self, seq: str) -> float:
        seq = seq.upper()
        if not seq:
            return 0.0
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq) * 100

    def find_simple_repeats(self, seq: str, max_len: int = 6) -> list:
        """查找长度在2~max_len之间的重复子串（简单实现）"""
        repeats = set()
        seq_upper = seq.upper()
        for length in range(2, max_len + 1):
            for i in range(len(seq_upper) - length + 1):
                sub = seq_upper[i:i+length]
                if seq_upper.count(sub) > 1:
                    repeats.add(sub)
        return sorted(list(repeats), key=lambda x: len(x), reverse=True)

    def calculate_tm(self, seq: str) -> float:
        """简单Tm公式：Tm = 4*(G+C) + 2*(A+T)（适用于短于14bp）"""
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        at = len(seq) - gc
        return 4 * gc + 2 * at

    def find_restriction_sites(self, seq: str) -> list:
        """查找常见酶切位点，返回 (酶名, 识别序列, 起始位置) 列表"""
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
