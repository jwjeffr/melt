import ovito
import sys
from dataclasses import dataclass
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt


@dataclass
class ProportionCalculator:

    """
    class for calculating proportion of structure (FCC, BCC, etc.) at a given frame of the simulation
    """

    structure: str
    num_atoms: int
    
    def __call__(self, frame: ovito.data.DataCollection) -> float:
    
        for key, value in frame.attributes.items():
            if key != f'PolyhedralTemplateMatching.counts.{self.structure}':
                continue
            proportion_solid = value / self.num_atoms
            
        return proportion_solid


def get_solid_percentage(input_file_name: str, cutoff: float, structure='FCC') -> float:

    """
    function for getting solid percentage at the final frame
    """
    
    pipeline = ovito.io.import_file(input_file_name)
    pipeline.modifiers.append(ovito.modifiers.PolyhedralTemplateMatchingModifier(rmsd_cutoff=cutoff))
    num_atoms = 0
    initial_frame = pipeline.compute(0)
    for key, value in initial_frame.attributes.items():
        if 'PolyhedralTemplateMatching.counts' not in key:
            continue
        num_atoms += value
    proportion_calculator = ProportionCalculator(structure=structure, num_atoms=num_atoms)
    
    final_frame = pipeline.compute(pipeline.source.num_frames)
    return proportion_calculator(final_frame)
        
if __name__ == '__main__':

    # take command line args
    # plot is outputted for each pressure
    # so if i wanted the percent solid-temperature curve for P = 1 bar
    # for temperatures between 1000 and 2000 K, with a step of 50 K, i would run:
    # python percent_solid.py 1000 2000 50 1

    low_temp = int(sys.argv[1])
    high_temp = int(sys.argv[2])
    step_temp = int(sys.argv[3])
    pressure = int(sys.argv[4])
    
    temps = np.arange(low_temp, high_temp + step_temp, step=step_temp)
    percent_solid = np.zeros(temps.shape)
    
    for index, temp in enumerate(temps):
    
        file_name = f'equil_{temp:.0f}_{pressure:.0f}.dump'
        solid_percentage = get_solid_percentage(file_name, cutoff=0.15)
        percent_solid[index] = solid_percentage
        print('done')
        
    plt.scatter(temps, percent_solid * 100, facecolor='lightblue', edgecolor='black', zorder=6)
    plt.xlabel('temperature (K)')
    plt.ylabel('final percent fcc after 0.5 ns')
    plt.grid()
    plt.savefig(f'melting_{pressure:.0f}.png', dpi=800, bbox_inches='tight')
