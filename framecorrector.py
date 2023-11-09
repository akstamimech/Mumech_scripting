#frame corrector
import os
from pathlib import Path

# Set the directory path
directory = Path.cwd() / 'frames_centered'

# Get a list of all .pdb files in the directory
files = [file for file in os.listdir(directory) if file.endswith('.pdb')]

# Iterate over the files
for file in files:
    # Extract the frame number from the file name
    frame_number = file.split('frame')[1].split('.pdb')[0]
    
    # Check if the frame number is less than 10
    if len(frame_number) == 2 and frame_number[0] == '0':
        # Create the new file name with single-digit frame number
        new_file_name = file.replace(frame_number, frame_number[1])
    else:
        # Convert the frame number to an integer
        frame_number = int(frame_number)
        
        # Check if the frame number is greater than or equal to 9000
        if frame_number >= 9000:
            # Convert the frame number to the new index
            new_frame_number = frame_number - 9000 + 90
            
            # Create the new file name with updated frame number
            new_file_name = file.replace(f'frame{frame_number}', f'frame{new_frame_number}')
        else:
            # Skip files with frame numbers between 10 and 89
            continue
    
    # Create the full paths for the old and new file names
    old_file_path = os.path.join(directory, file)
    new_file_path = os.path.join(directory, new_file_name)
    
    # Rename the file
    os.rename(old_file_path, new_file_path)