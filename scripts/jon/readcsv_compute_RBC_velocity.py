import sys
sys.path.append("L:/hemocell/scripts/measure/")  # 'default'
import HCELL_readcsv
import numpy as np

DX = 5e-7;   DT = 1e-7

BEGIN =   100000   # 1e5 / one-hundred thousand
END  =  10000000   # 1e7 / ten million  (inclusive)
STEP =    100000
#%%

# RBCs can be still handled by preinlet but their center is past the preinlet border
# Subtract this slack variable to the preinlet border (410.5) to include these
SLACK = 7.5  

# Allow for allow measuring some RBCs / Platelets depending on their x-location
CONSTRAINT1 = lambda x : True if x < 410.5 - SLACK else False   # Entire geometry without preinlet
CONSTRAINT2 = lambda x : True if x > 410.5 - SLACK else False   # Preinlet
CONSTRAINT3 = lambda x : True if x < 410.5 - SLACK and x > 363.5 else False   # Inlet (excluding preinlet)
CONSTRAINT4 = lambda x : True if x < 36.5  else False              # Outlet

CONSTRAINTS = [CONSTRAINT1, CONSTRAINT2, CONSTRAINT3, CONSTRAINT4]

cNames = ["geometry", "preinlet", "inlet", "outlet"]

#%%
for datapath in ["F:/buggy/AR1/output/csv/",       "F:/buggy/AR2/output/csv/", \
                 "F:/buggy/AR2_stiff/output/csv/", "F:/buggy/AR3/output/csv/", ]:
    
    for CONSTRAINT, cname in zip(CONSTRAINTS, cNames):
        print(datapath + " " + cname)

        # Read RBC data
        RBC, PLT, _ = HCELL_readcsv.open_csv_files(r=True, p=True, ct3=False, begin=BEGIN, \
                                                   end=END+STEP, timestep=STEP, datapath=datapath, \
                                                   rbcname="RBC", pltname="PLT")
        
        # Compute velocities for RBCs and platelets
        for CELL, label in zip([RBC, PLT], ["RBC", "PLT"]):
            
            total_cell_count_over_time = 0
            T = len(CELL)
            avg_vel_magnitudes_over_time = []
            std_vel_magnitudes_over_time = []
            
            # Iterate over time
            for t, cell_t_instance in enumerate(CELL):                
                # Put cell IDs and velocities into numpy array once
                positions  = np.array(cell_t_instance.position)
                velocities = np.array(cell_t_instance.velocity)
                
                # Compute average velocity over cells
                vel_magnitudes = []
                cell_count = 0
                # Iterate over cells
                for pos_triplet, vel_triplet in zip(positions, velocities):
                    # Only look at cells in inlet / outlet
                    if CONSTRAINT(pos_triplet[0]):            
                        cell_count += 1
                        # Compute velocity magnitude from its x,y,z components
                        vel_magnitude = np.linalg.norm(vel_triplet)
                        vel_magnitudes.append(vel_magnitude)
      
                if cell_count != 0:
                    # Average velocity for this timestep
                    avg_vel_magnitude = np.mean(vel_magnitudes)
                    avg_vel_magnitudes_over_time.append(avg_vel_magnitude)
                    
                     # stddev of velocities over this timestep
                    std_vel_magnitude = np.std(vel_magnitudes)
                    std_vel_magnitudes_over_time.append(np.std(vel_magnitudes))
                    
                    total_cell_count_over_time += cell_count
            
            # Average velocity over time, over all cells per time instance
            avg_avg_vel_magnitude = np.mean(avg_vel_magnitudes_over_time)
            # Standard deviation over time
            #std_avg_vel_magnitude = np.std( avg_vel_magnitudes_over_time)
                
            # Average over time, convert to mm/s, print
            conversion_factor = DX / DT * 1000   #  * 1000 for m/s -> mm/s
            
            # Make final statistics
            avg_speed = avg_avg_vel_magnitude * conversion_factor
            avg_sd_over_cells = np.mean(std_vel_magnitudes_over_time) * conversion_factor
            sd_over_time = np.std(avg_vel_magnitudes_over_time) * conversion_factor
            
            print(label)
            #print(f"Number of cells iterated over: {total_cell_count_over_time}")
            #print(f"Average speed: {avg_speed:.2f}")
            #print(f"Average std dev over cells in one time instance: {avg_sd_over_cells:.2f}")
            #print(f"std dev of average velocities over time: {sd_over_time:.2f}")
            print(f"{avg_speed:.2f} {total_cell_count_over_time} {avg_sd_over_cells:.2f} {sd_over_time:.2f}")
            print()