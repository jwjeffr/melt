import ovito
import sys
from dataclasses import dataclass
import numpy as np
from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt


@dataclass
class ProportionCalculator:

    structure: str
    num_atoms: int
    
    def __call__(self, frame: ovito.data.DataCollection) -> float:
    
        for key, value in frame.attributes.items():
            if key != f'PolyhedralTemplateMatching.counts.{self.structure}':
                continue
            proportion_solid = value / self.num_atoms
            
        return proportion_solid


def get_solid_percentage(input_file_name: str, structure='FCC', cutoff=None) -> float:
    
    tolerance = 0.01
    proportion_solid = 1.0
    proportion_calculator = None
    
    if cutoff is None:
    
        cutoff = 0.15
    
        # need to select an appropriate cutoff such that the solid-(solid + liquid) proportion is about 55% initially
        
        while proportion_solid > 0.55 - tolerance:
        
            # import file, append polyhedral template matching for structure identification
        
            pipeline = ovito.io.import_file(input_file_name)
            pipeline.modifiers.append(ovito.modifiers.PolyhedralTemplateMatchingModifier(rmsd_cutoff=cutoff))
            
            initial_frame = pipeline.compute(0)
            
            if not proportion_calculator:
                num_atoms = 0
                for key, value in initial_frame.attributes.items():
                    if 'PolyhedralTemplateMatching.counts' not in key:
                        continue
                    num_atoms += value
                proportion_calculator = ProportionCalculator(structure=structure, num_atoms=num_atoms)
                    
            proportion_solid = proportion_calculator(initial_frame)
            cutoff += -0.01
            
    else:
    
        pipeline = ovito.io.import_file(input_file_name)
        pipeline.modifiers.append(ovito.modifiers.PolyhedralTemplateMatchingModifier(rmsd_cutoff=cutoff))
        num_atoms = 0
        initial_frame = pipeline.compute(0)
        for key, value in initial_frame.attributes.items():
            if 'PolyhedralTemplateMatching.counts' not in key:
                continue
            num_atoms += value
        proportion_calculator = ProportionCalculator(structure=structure, num_atoms=num_atoms)
            
    # calculate final solid-(solid + liquid) proportion
    
    final_frame = pipeline.compute(pipeline.source.num_frames)
    return proportion_calculator(final_frame)
    
    
def logistic(temp, normalization, center):

    return 1.0 / (1.0 + np.exp(-(temp - center) / normalization))
        
if __name__ == '__main__':

    low_temp = float(sys.argv[1])
    high_temp = float(sys.argv[2])
    step_temp = float(sys.argv[3])
    
    temps = np.arange(low_temp, high_temp + step_temp, step=step_temp)
    percent_solid = np.zeros(temps.shape)
    pressure = 1
    
    for index, temp in enumerate(temps):
    
        file_name = f'equil_{temp:.0f}_{pressure:.0f}.dump'
        solid_percentage = get_solid_percentage(file_name, cutoff=0.15)
        percent_solid[index] = solid_percentage
        print('done')
        
    plt.scatter(temps, percent_solid * 100, facecolor='lightblue', edgecolor='black', zorder=6)
    plt.xlabel('temperature (K)')
    plt.ylabel('final percent fcc after 0.5 ns')
    plt.grid()
    plt.savefig('melting.png', dpi=800, bbox_inches='tight')
