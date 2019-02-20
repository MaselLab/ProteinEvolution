import numpy as np


'''
Jason's hydrophobic AA run lengths functions
'''
			
hydrophobicity_map={'FILVMW':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':1, 
                    'A':-1,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':-1,'H':-1,
                    'K':-1,'P':-1,'S':-1,'T':-1,'Y':-1,'X':0,'U':-1},
                    'FILVM':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':-1, 
                    'A':-1,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':-1,'H':-1,
                    'K':-1,'P':-1,'S':-1,'T':-1,'Y':-1,'X':0,'U':-1},
                    'AGWPY':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':0, 
                    'A':0,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':0,'H':-1,
                    'K':-1,'P':0,'S':-1,'T':-1,'Y':0,'X':0,'U':-1},
                    'FLIMVWAG':{'F':1,'I':1,'L':1,'V':1,'M':1,'W':1, 
                    'A':1,'R':-1,'N':-1,'D':-1,'C':-1,'E':-1,'Q':-1,'G':1,'H':-1,
                    'K':-1,'P':-1,'S':-1,'T':-1,'Y':-1,'X':0,'U':-1}}
			
def run_lengths_phobic(sequence_hydro):
    counts=np.zeros(len(sequence_hydro))
    first_hydrophobic=next(x for x,y in enumerate(sequence_hydro) if y==1)
    prev=1
    running_count=0
    phobicity_flag=1
    for amino in sequence_hydro[first_hydrophobic+1:]:
        if amino != prev:
            if phobicity_flag>0:
                counts[running_count]=counts[running_count]+1
                
            phobicity_flag=-1*phobicity_flag
            running_count=0
        else:
            running_count=running_count+1
            
        prev=amino
    
    if phobicity_flag>0:
                counts[running_count]=counts[running_count]+1    
    
    return counts  
	
def run_lengths_phobic_metric(sequence_amino,hydro_map):
    sequence_hydro=[hydrophobicity_map[hydro_map][_] for _ in sequence_amino]
    total_hydro=sum([_ for _ in sequence_hydro if _==1])
    if total_hydro==0:
        return 0.
    counts=run_lengths_phobic(sequence_hydro)
    meanRL=total_hydro/sum(counts)
    hp=float(total_hydro)/len(sequence_hydro)
    norm=1/(1-hp)
    return meanRL/norm
