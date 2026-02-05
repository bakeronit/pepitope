from dataclasses import dataclass

@dataclass
class VaccineType:
    subtype: str
    name: str
    sequence: str
    ha1_region: tuple[int, int]
    epitope_positions: dict[list[int]]
    ve_equation_params: dict[dict[float]]

    def __str__(self):
        return f"{self.name}\n{self.sequence}\nHA1 protein region:{self.ha1_region}\nEpitope positions:\n{self.epitope_positions}\nVE equation parameters:\n{self.ve_equation_params}"

    def __repr__(self):
        return self.__str__()

H3N2 = VaccineType(
    subtype = "H3N2",
    name = "A/Hong Kong/4801/2014 (EPI614406)",
    sequence = "MKTIIALSYILCLVFAQKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGEICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCKSECITPNGSIPNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGRGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHNVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI",
    ha1_region = (17, 344),
    epitope_positions = {
        "A" : [122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168],
        "B" : [128, 129, 155, 156, 157, 158, 159, 160, 163, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198],
        "C" : [44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312],
        "D" : [96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248],
        "E" : [57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265]
    },
    ve_equation_params = {
        ("CDC", "Efficacy"     ) : { "m": -3.0580, "m_err": 0.50132, "b": 0.59394, "b_err": 0.072578 },
        ("CDC", "Effectiveness") : { "m": -3.3218, "m_err": 0.55778, "b": 0.65505, "b_err": 0.080753 },
        ("NH", "Efficacy"      ) : { "m": -2.6132, "m_err": 0.38326, "b": 0.51604, "b_err": 0.049732 },
        ("NH", "Effectiveness" ) : { "m": -2.6614, "m_err": 0.53525, "b": 0.57682, "b_err": 0.069455 }
    }
)