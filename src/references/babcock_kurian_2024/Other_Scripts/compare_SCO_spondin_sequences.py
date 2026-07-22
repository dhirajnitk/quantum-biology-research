import numpy as np
from difflib import SequenceMatcher

MAPPING = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     	   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     	   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     	   'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def shorten(three_letter_code):
	return MAPPING[three_letter_code]

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()
    
def similar_substring(smaller_str, bigger_str):
	similarities = []
	for i in range(len(bigger_str) - 2*len(smaller_str)):
		similarities.append(similar(smaller_str, bigger_str[i:len(smaller_str)+i]))

	similarities = np.array(similarities)
	indmax = np.argmax(similarities)
	return [indmax, similarities, bigger_str[indmax:indmax+len(smaller_str)]]


with open("PDB_files/P98167.fasta", 'r') as f:
	whole_sequence = f.readlines()[1:]
	whole_sequence = [line[:-1] for line in whole_sequence]
	whole_sequence = "".join(whole_sequence)

	
with open("PDB_files/AF-O02660-F1-model_v4.pdb", 'r') as f:
	model_002660 = f.readlines()
	model_002660 = [line.split()[4:] for line in model_002660 if line[:6] == "SEQRES"]
	model_002660 = [shorten(aa) for arr in model_002660 for aa in arr]
	model_002660 = "".join(model_002660)

with open("PDB_files/AF-O02661-F1-model_v4.pdb", 'r') as f:
	model_002661 = f.readlines()
	model_002661 = [line.split()[4:] for line in model_002661 if line[:6] == "SEQRES"]
	model_002661 = [shorten(aa) for arr in model_002661 for aa in arr]
	model_002661 = "".join(model_002661)

with open("PDB_files/AF-Q8SPM5-F1-model_v4.pdb", 'r') as f:
	model_q8spm5 = f.readlines()
	model_q8spm5 = [line.split()[4:] for line in model_q8spm5 if line[:6] == "SEQRES"]
	model_q8spm5 = [shorten(aa) for arr in model_q8spm5 for aa in arr]
	model_q8spm5 = "".join(model_q8spm5)

with open("PDB_files/AF-Q9TTS5-F1-model_v4.pdb", 'r') as f:
	model_q9tts5 = f.readlines()
	model_q9tts5 = [line.split()[4:] for line in model_q9tts5 if line[:6] == "SEQRES"]
	model_q9tts5 = [shorten(aa) for arr in model_q9tts5 for aa in arr]
	model_q9tts5 = "".join(model_q9tts5)

with open("PDB_files/AF-Q9TTS4-F1-model_v4.pdb", 'r') as f:
	model_q9tts4 = f.readlines()
	model_q9tts4 = [line.split()[4:] for line in model_q9tts4 if line[:6] == "SEQRES"]
	model_q9tts4 = [shorten(aa) for arr in model_q9tts4 for aa in arr]
	model_q9tts4 = "".join(model_q9tts4)

alphafold_models = [model_q8spm5, model_q9tts4, model_q9tts5, model_002661, model_002660]
mdls_params = []

seq = "TGVCTAGCACPTGLFLHNSSCLPPSQCPCQLRGQLYAPGAVARLDSCSNCTCISGEMVCASEPCPVACGWSPWTPWSLCSRSCNVGVRRRFRAGTAPPAAFGGAACQGPNMEAEFCSLRPCGGPAGEWGPWSPCSVPCGGGYRNRTRGSSGPSPVDFSTCGLQPCAGPAPGVCPPGKRWLDCAQGPASCAELSAPRGADQPCHPGCYCPSGMLLLNNACVPTQDCPCTHGGRLHPPGSAVLRPCENCSCVSGLITNCTSWPCKEGQPTWSPWTPWSECSASCGPARRHKHRFCTRPPGGAPSSMAPPLLLSSVPPLCPGPEAEEEPCLLPECDRAGGWGPWGPWSSCSRSCGGGLRSRSRACDQPPPQGLGDYCEGPRAQGAACQALPCPVTNCTAIEGAEYSACGPPCPRSCDDLVHCVWHCQPGCYCPPGQVLSADGTVHVQPGHCSCLDLLTGERHRPGAQLAKPDGCNYCTCSEGQLTCTDLPCPVPGAWCPWSEWTACSQPCQGQTRTRSRACSCPAPQHGGAPCPGEAGEAGAQHQRETCASTPECPVDGAWSPWGPWSPCEVCLGRSHRSRECSWPPTSEGGRPCPGGHRQSRPCQGNSTQCTDCAGGQDLLPCGQPCPRSCEDLSPGVECQPDSMGCQQPRCGCPPGQLSQDGLCVTPSQCRCQYQPGAMGIPENQS"

'''
for mdl in alphafold_models:
	indmax, similarities, _ = similar_substring(mdl, whole_sequence)
	print(mdl)
	print(similarities[indmax])
	print()
'''

cnt = 0
for c1, c2 in zip(seq, alphafold_models[2]):
	if c1 != c2:
		cnt = cnt + 1
		print("___", end="")
	else:
		print(c1, end="")

print()
print(cnt)
#print(len(seq))
#print(len(alphafold_models[2]))

'''
0 1.0
346 1.0
2897 0.9985401459854014
2561 1.0
889 1.0

'''
#print(similarities[indmax])
#print(f"{whole_sequence[0:indmax]}__{whole_sequence[indmax:indmax+len(model_002660)]}__{whole_sequence[indmax+len(model_002660):]}")

#print("f[__{whole_sequence[0:len(model_q8spm5)]}__]{whole_sequence[len(model_q8spm5):346])[__{whole_sequence[346:len(model_q9tts4)+346]}__]{whole_sequence

'''
with open("AF-O02661-F1-model_v4.pdb", 'r') as f:
with open("AF-Q8SPM5-F1-model_v4.pdb", 'r') as f:
with open("AF-Q9TTS4-F1-model_v4.pdb", 'r') as f:
with open("AF-Q9TTS5-F1-model_v4.pdb", 'r') as f:
'''
