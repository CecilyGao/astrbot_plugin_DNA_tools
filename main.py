import re
import math
from typing import List, Optional, Tuple
from astrbot.api.event import filter, AstrMessageEvent
from astrbot.api.star import Context, Star, register
from astrbot.api import logger, AstrBotConfig

# ======================== BioPython 导入 ========================
try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC, MeltingTemp as mt
    from Bio.Restriction import Restriction
    from Bio.Data import CodonTable
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logger.warning("BioPython 未安装，部分功能将使用简化算法。")

# ======================== 常量 ========================
ALLOWED_BASES = set("ATCGUN")

# ======================== 常用限制酶列表（名称 -> 识别序列） ========================
# 包含 4bp、6bp、8bp 识别位点的常见酶，供回退模式使用
# 注：BioPython 可用时，将直接使用其酶对象，该字典仅作为备用字符串查找
COMMON_ENZYMES = {
    # 4bp 切点酶
    'AluI': 'AGCT',
    'HaeIII': 'GGCC',
    'MspI': 'CCGG',
    'TaqI': 'TCGA',
    'RsaI': 'GTAC',
    'HinfI': 'GANTC',      # 含简并
    'MboI': 'GATC',
    'Sau3AI': 'GATC',
    'DpnI': 'GATC',        # 甲基化敏感，但此处仅识别序列
    'BstUI': 'CGCG',
    'HpaII': 'CCGG',
    'MaeI': 'CTAG',
    'MnlI': 'CCTC',
    # 6bp 切点酶
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'NotI': 'GCGGCCGC',    # 8bp
    'XbaI': 'TCTAGA',
    'SalI': 'GTCGAC',
    'PstI': 'CTGCAG',
    'SmaI': 'CCCGGG',
    'KpnI': 'GGTACC',
    'SacI': 'GAGCTC',
    'NcoI': 'CCATGG',
    'NdeI': 'CATATG',
    'XhoI': 'CTCGAG',
    'SpeI': 'ACTAGT',
    'BglII': 'AGATCT',
    'MfeI': 'CAATTG',
    'ClaI': 'ATCGAT',
    'HpaI': 'GTTAAC',
    'MluI': 'ACGCGT',
    'NheI': 'GCTAGC',
    'AgeI': 'ACCGGT',
    'BsrGI': 'TGTACA',
    'BstEII': 'GGTNACC',   # 含简并
    'BstXI': 'CCANNNNNNTGG',
    'DraI': 'TTTAAA',
    'EcoRV': 'GATATC',
    'HincII': 'GTYRAC',    # 含简并
    'PvuII': 'CAGCTG',
    'ScaI': 'AGTACT',
    'SphI': 'GCATGC',
    'StuI': 'AGGCCT',
    'AatII': 'GACGTC',
    'BglI': 'GCCNNNNNGGC',
    'NruI': 'TCGCGA',
    'SfiI': 'GGCCNNNNNGGCC', # 8bp
    'BspHI': 'TCATGA',
    'XcmI': 'CCANNNNNNNNNTGG',
    # 更多...
}

# ======================== 主插件类 ========================
@register(
    "astrbot_plugin_DNA_tools",
    "CecilyGao",
    "参考擎科生物官网云工具，使用 BioPython 实现更准确的 DNA 分析",
    "0.0.1",
    "https://github.com/CecilyGao/astrbot_plugin_DNA_tools"
)
class DNAPlugin(Star):
    def __init__(self, context: Context, config: AstrBotConfig):
        super().__init__(context)
        self.config = config
        self.group_whitelist = config.get('group_whitelist', [])
        self.format_width = config.get('format_width', 60)
        # Tm 计算参数（可配置）
        self.salt_conc = config.get('salt_conc', 0.05)      # 单价阳离子浓度 (M)
        self.primer_conc = config.get('primer_conc', 0.2)   # 引物浓度 (μM)

    # ======================== 权限辅助 ========================
    def _check_group_permission(self, event: AstrMessageEvent) -> bool:
        if not self.group_whitelist:
            return True
        group_id = str(event.get_group_id())
        return group_id in self.group_whitelist

    # ======================== 主命令入口 ========================
    @filter.command("DNA", alias={'dna'})
    async def dna_command(self, event: AstrMessageEvent):
        if not self._check_group_permission(event):
            yield event.plain_result("❌ 当前群未在 DNA 工具白名单中。")
            return

        full_text = event.message_str.strip()
        parts = full_text.split()
        if len(parts) < 2:
            yield event.plain_result(self._get_help())
            return

        subcmd = parts[1].lower()
        args = parts[2:] if len(parts) > 2 else []

        if subcmd in ['help', '帮助']:
            yield event.plain_result(self._get_help())
            return

        seq_required_cmds = ['format', '格式化', 'reverse', '反向互补', 
                             'translate', '翻译', 'gc', 'gc含量',
                             'primer', '引物分析', 'restriction', '酶切']
        if subcmd in seq_required_cmds:
            seq = self._extract_sequence(args)
            if seq is None:
                yield event.plain_result("⚠️ 请提供有效的 DNA 序列（只允许 A、T、C、G、U、N），例如：`DNA format ATCG`")
                return
        else:
            yield event.plain_result(f"❌ 未知子命令：{subcmd}\n" + self._get_help())
            return

        # 子命令分发
        if subcmd in ['format', '格式化']:
            result = await self._cmd_format(seq)
        elif subcmd in ['reverse', '反向互补']:
            result = await self._cmd_reverse(seq)
        elif subcmd in ['translate', '翻译']:
            result = await self._cmd_translate(seq)
        elif subcmd in ['gc', 'gc含量']:
            result = await self._cmd_gc(seq)
        elif subcmd in ['primer', '引物分析']:
            result = await self._cmd_primer(seq)
        elif subcmd in ['restriction', '酶切']:
            result = await self._cmd_restriction(seq)
        else:
            result = f"❌ 未知子命令：{subcmd}\n" + self._get_help()
        yield event.plain_result(result)

    # ======================== 参数提取 ========================
    @staticmethod
    def _extract_sequence(args: List[str]) -> Optional[str]:
        if not args:
            return None
        raw = ''.join(args).replace(' ', '').replace('\n', '').replace('\r', '')
        seq = raw.upper()
        if not seq:
            return None
        if not set(seq).issubset(ALLOWED_BASES):
            return None
        return seq

    # ======================== 各子命令实现 ========================
    async def _cmd_format(self, seq: str) -> str:
        width = self.format_width
        formatted = '\n'.join([seq[i:i+width] for i in range(0, len(seq), width)])
        return f"📄 格式化结果（每行 {width} 个碱基）：\n{formatted}"

    async def _cmd_reverse(self, seq: str) -> str:
        if BIOPYTHON_AVAILABLE:
            try:
                rev_comp = str(Seq(seq).reverse_complement())
            except Exception as e:
                logger.error(f"BioPython reverse complement error: {e}")
                rev_comp = self._fallback_reverse_complement(seq)
        else:
            rev_comp = self._fallback_reverse_complement(seq)
        return f"🧬 反向互补序列：\n{rev_comp}"

    async def _cmd_translate(self, seq: str) -> str:
        if BIOPYTHON_AVAILABLE:
            try:
                codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
                protein = str(Seq(seq).translate(table=codon_table, stop_symbol="*"))
            except Exception as e:
                logger.error(f"BioPython translate error: {e}")
                protein = self._fallback_translate(seq)
        else:
            protein = self._fallback_translate(seq)
        return f"🧪 翻译结果（标准密码子表）：\n{protein}"

    async def _cmd_gc(self, seq: str) -> str:
        if BIOPYTHON_AVAILABLE:
            try:
                gc = GC(seq)
            except Exception as e:
                logger.error(f"BioPython GC error: {e}")
                gc = self._fallback_gc(seq)
        else:
            gc = self._fallback_gc(seq)
        repeats = self._find_simple_repeats(seq)
        repeat_str = ', '.join(repeats[:5]) if repeats else "无"
        return f"📊 GC 含量：{gc:.2f}%\n🔁 重复序列（长度2-6，仅显示前5个）：{repeat_str}"

    async def _cmd_primer(self, seq: str) -> str:
        # 使用 BioPython 计算 Tm（最近邻法，含盐校正）
        if BIOPYTHON_AVAILABLE:
            try:
                tm = mt.Tm_NN(Seq(seq),
                              Na=self.salt_conc,
                              K=0.0,
                              Tris=0.0,
                              Mg=0.0,
                              dNTPs=0.0,
                              salt_correction=mt.salt_correction_schildkraut,
                              DNA=True,
                              self_comp=False)
            except Exception as e:
                logger.error(f"BioPython Tm calculation error: {e}")
                tm = self._fallback_tm(seq)
        else:
            tm = self._fallback_tm(seq)

        gc = self._gc_content(seq)  # 统一用我们自己的 GC 计算（也可直接用 BioPython）
        return (f"🧬 引物分析（使用 BioPython 最近邻法，盐浓度 {self.salt_conc} M，引物浓度 {self.primer_conc} μM）：\n"
                f"序列：{seq}\n"
                f"GC 含量：{gc:.2f}%\n"
                f"Tm 值（修正后）：{tm:.2f} °C")

    async def _cmd_restriction(self, seq: str) -> str:
        """
        查找限制酶切位点。
        若 BioPython 可用，则使用其酶对象进行模糊匹配（支持简并碱基）；
        否则回退到内置字典进行精确字符串查找。
        """
        if BIOPYTHON_AVAILABLE:
            try:
                seq_obj = Seq(seq)
                found = []
                # 遍历常用酶列表
                for enzyme_name in COMMON_ENZYMES.keys():
                    try:
                        enzyme = getattr(Restriction, enzyme_name)
                    except AttributeError:
                        # 如果 BioPython 中没有该酶，则跳过
                        continue
                    # 搜索该酶的所有切割位点（返回位置列表）
                    sites = enzyme.search(seq_obj)
                    if sites:
                        # 获取识别序列（可能含简并碱基）
                        site_seq = str(enzyme.site)
                        for pos in sites:
                            found.append((enzyme_name, site_seq, pos))
                if found:
                    lines = ["🔬 发现的酶切位点（使用 BioPython 搜索）："]
                    for name, site, pos in found:
                        lines.append(f"  {name} 识别序列 {site} 切割位置 {pos}")
                    return "\n".join(lines)
                else:
                    return "❌ 未发现已知酶切位点。"
            except Exception as e:
                logger.error(f"BioPython restriction error: {e}")
                return self._fallback_restriction(seq)
        else:
            return self._fallback_restriction(seq)

    # ======================== 回退函数（BioPython 不可用时使用） ========================
    def _gc_content(self, seq: str) -> float:
        """统一 GC 含量计算（优先使用 BioPython）"""
        if BIOPYTHON_AVAILABLE:
            try:
                return GC(seq)
            except:
                return self._fallback_gc(seq)
        else:
            return self._fallback_gc(seq)

    def _fallback_gc(self, seq: str) -> float:
        seq = seq.upper()
        if not seq:
            return 0.0
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq) * 100

    def _fallback_reverse_complement(self, seq: str) -> str:
        complement = {'A':'T','T':'A','C':'G','G':'C','U':'A','N':'N'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def _fallback_translate(self, seq: str) -> str:
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
        return ''.join(protein)

    def _fallback_tm(self, seq: str) -> float:
        # 简单 Wallace 公式（短引物适用）
        gc = seq.upper().count('G') + seq.upper().count('C')
        at = len(seq) - gc
        return 4 * gc + 2 * at

    def _fallback_restriction(self, seq: str) -> str:
        """回退模式：使用内置字典进行精确字符串查找（不支持简并碱基）"""
        found = []
        seq_upper = seq.upper()
        for enzyme, site in COMMON_ENZYMES.items():
            pos = seq_upper.find(site)
            if pos != -1:
                found.append((enzyme, site, pos))
        if found:
            lines = ["🔬 发现的酶切位点（使用内置字典精确匹配）："]
            for enzyme, site, pos in found:
                lines.append(f"  {enzyme} 识别序列 {site} 位置 {pos}")
            return "\n".join(lines)
        else:
            return "❌ 未发现已知酶切位点。"

    def _find_simple_repeats(self, seq: str, max_len: int = 6) -> List[str]:
        repeats = set()
        seq_upper = seq.upper()
        for length in range(2, max_len + 1):
            for i in range(len(seq_upper) - length + 1):
                sub = seq_upper[i:i+length]
                if seq_upper.count(sub) > 1:
                    repeats.add(sub)
        return sorted(list(repeats), key=lambda x: len(x), reverse=True)

    # ======================== 帮助信息 ========================
    def _get_help(self) -> str:
        help_text = (
            "🧬 DNA 工具使用帮助（基于 BioPython 增强版）\n"
            "命令格式：DNA <子命令> [序列]\n\n"
            "可用子命令：\n"
            "  format 或 格式化   <序列>  – 每行60碱基格式化输出\n"
            "  reverse 或 反向互补 <序列>  – 计算反向互补序列（BioPython 实现）\n"
            "  translate 或 翻译  <序列>  – 标准密码子表翻译（BioPython 实现）\n"
            "  gc 或 gc含量       <序列>  – GC 含量（BioPython） + 简单重复序列\n"
            "  primer 或 引物分析 <序列>  – GC% + Tm 值（最近邻法，含盐校正，BioPython）\n"
            "  restriction 或 酶切 <序列> – 查找常用限制酶切位点（支持简并，BioPython）\n"
            "  help 或 帮助              – 显示此帮助\n\n"
            "示例：\n"
            "  DNA format ATCGATCGATCG\n"
            "  DNA 反向互补 ATGC\n"
            "  DNA translate ATGCCGTAA\n"
            "  DNA gc ATGCGC\n"
            "  DNA primer ATGCGT\n"
            "  DNA 酶切 GAATTC\n\n"
            "配置项（在插件配置文件中可调整）：\n"
            "  salt_conc: 单价阳离子浓度 (M)，默认 0.05\n"
            "  primer_conc: 引物浓度 (μM)，默认 0.2\n"
            "  format_width: 格式化每行碱基数，默认 60\n"
            "  group_whitelist: 允许使用的群组 ID 列表，留空则全部可用\n\n"
            "注意：序列中只允许 A、T、C、G、U、N，空格和换行自动忽略。\n" 
            "都用本插件了，想必也知道序列应该怎么输入吧。"
        )
        if not BIOPYTHON_AVAILABLE:
            help_text += "\n⚠️ BioPython 未安装，已回退至简化算法，建议安装以获得更准确结果。"
        return help_text