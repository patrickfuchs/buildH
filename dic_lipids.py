# Dicts for reconstructing hydrogens in POPC
# Format:
# dic= { "atom1": ("typeofH2build", "helper1", "helper2"),
#        "atom2": ("typeofH2build", "helper1", "helper2"),
#        ...}
# NOTE If you add a new lipid, we advice the user to add all the possible carbons on which we want to rebuild hydrogens. If the -opx option of buildH is used, it is strictly mandatory.

import copy

Berger_POPC = {
        # residue name
        "resname": "POPC",
        # choline
        "C1": ("CH3", "N4", "C5"),
        "C2": ("CH3", "N4", "C5"),
        "C3": ("CH3", "N4", "C5"),
        "C5": ("CH2", "N4", "C6"),
        "C6": ("CH2", "C5", "O7"),
        # glycerol
        "C12": ("CH2", "O11", "C13"),
        "C13": ("CH", "C12", "C32", "O14"),
        "C32": ("CH2", "C13", "O33"),
        # sn2
        "C17": ("CH2", "C15", "C18"),
        "C18": ("CH2", "C17", "C19"),
        "C19": ("CH2", "C18", "C20"),
        "C20": ("CH2", "C19", "C21"),
        "C21": ("CH2", "C20", "C22"),
        "C22": ("CH2", "C21", "C23"),
        "C23": ("CH2", "C22", "C24"),
         # C24=C25 --> double bond
        "C24": ("CHdoublebond", "C23", "C25"),
        "C25": ("CHdoublebond", "C24", "C26"),
        "C26": ("CH2", "C25", "C27"),
        "C27": ("CH2", "C26", "C28"),
        "C28": ("CH2", "C27", "C29"),
        "C29": ("CH2", "C28", "C30"),
        "C30": ("CH2", "C29", "C31"),
        "C31": ("CH2", "C30", "CA1"),
        "CA1": ("CH2", "C31", "CA2"),
        "CA2": ("CH3", "CA1", "C31"), # helper1 is the first C connected to CH3, helper is 2 atoms away
        # sn1
        "C36": ("CH2", "C34", "C37"),
        "C37": ("CH2", "C36", "C38"),
        "C38": ("CH2", "C37", "C39"),
        "C39": ("CH2", "C38", "C40"),
        "C40": ("CH2", "C39", "C41"),
        "C41": ("CH2", "C40", "C42"),
        "C42": ("CH2", "C41", "C43"),
        "C43": ("CH2", "C42", "C44"),
        "C44": ("CH2", "C43", "C45"),
        "C45": ("CH2", "C44", "C46"),
        "C46": ("CH2", "C45", "C47"),
        "C47": ("CH2", "C46", "C48"),
        "C48": ("CH2", "C47", "C49"),
        "C49": ("CH2", "C48", "C50"),
        "C50": ("CH3", "C49", "C48") # helper1 is the first C connected to CH3, helper is 2 atoms away
        }

# Alternative name for POPC in Berger lipids.
Berger_PLA = copy.deepcopy(Berger_POPC)
Berger_PLA["resname"] = "PLA"

Berger_POP = copy.deepcopy(Berger_POPC)
Berger_POP["resname"] = "POP"

# CHARMM POPC.
CHARMM_POPC = {
        # residue name
        "resname": "POPC",
        # choline
        "C13": ("CH3", "N", "C12"),
        "C14": ("CH3", "N", "C12"),
        "C15": ("CH3", "N", "C12"),
        "C12": ("CH2", "N", "C11"),
        "C11": ("CH2", "C12", "O12"),
        # glycerol
        "C1": ("CH2", "O11", "C2"),
        "C2": ("CH", "C1", "C3", "O21"),
        "C3": ("CH2", "C2", "O31"),
        # sn2
        "C22": ("CH2", "C21", "C23"),
        "C23": ("CH2", "C22", "C24"),
        "C24": ("CH2", "C23", "C25"),
        "C25": ("CH2", "C24", "C26"),
        "C26": ("CH2", "C25", "C27"),
        "C27": ("CH2", "C26", "C28"),
        "C28": ("CH2", "C27", "C29"),
         # C29=C210 --> double bond
        "C29": ("CHdoublebond", "C28", "C210"),
        "C210": ("CHdoublebond", "C29", "C211"),
        "C211": ("CH2", "C210", "C212"),
        "C212": ("CH2", "C211", "C213"),
        "C213": ("CH2", "C212", "C214"),
        "C214": ("CH2", "C213", "C215"),
        "C215": ("CH2", "C214", "C216"),
        "C216": ("CH2", "C215", "C217"),
        "C217": ("CH2", "C216", "C218"),
        "C218": ("CH3", "C217", "C216"), # helper1 is the first C connected to CH3, helper is 2 atoms away
        # sn1
        "C32": ("CH2", "C31", "C33"),
        "C33": ("CH2", "C32", "C34"),
        "C34": ("CH2", "C33", "C35"),
        "C35": ("CH2", "C34", "C36"),
        "C36": ("CH2", "C35", "C37"),
        "C37": ("CH2", "C36", "C38"),
        "C38": ("CH2", "C37", "C39"),
        "C39": ("CH2", "C38", "C310"),
        "C310": ("CH2", "C39", "C311"),
        "C311": ("CH2", "C310", "C312"),
        "C312": ("CH2", "C311", "C313"),
        "C313": ("CH2", "C312", "C314"),
        "C314": ("CH2", "C313", "C315"),
        "C315": ("CH2", "C314", "C316"),
        "C316": ("CH3", "C315", "C314") # helper1 is the first C connected to CH3, helper is 2 atoms away
        }
