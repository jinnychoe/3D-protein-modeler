# DNA Sequence Analysis and Protein 3D Modeling Tool

# The "DNA Sequence Analysis and Protein 3D Modeling Tool" is a Python program designed for molecular biology and bioinformatics enthusiasts. 
# This versatile tool enables users to perform various operations on DNA sequences, including counting nucleotides, converting DNA to RNA, and translating DNA into a protein sequence. 
# Additionally, it provides the functionality to generate 3D predictive models of protein structures based on DNA or user-inputted protein sequences.

import requests
from io import StringIO
import pymol
import os 
import gzip 
import subprocess 

# Create DNA sequence class
class DNA:
    def __init__(self, sequence): 
        self.sequence = sequence # Initialize variable

    def count_nucleotides(self): # Function counts nucleotides in DNA sequence object
        counts = [0, 0, 0, 0]  # Initializes nucleotide counts with list
        nucleotides = ['A', 'T', 'G', 'C']  # Defines order of nucleotides in list
        for nucleotide in self.sequence:
            if nucleotide in nucleotides:
                index = nucleotides.index(nucleotide)
                counts[index] += 1  # Accumulator tracks counts of each nucleotide
        return tuple(counts)  # A tuple of nucleotide counts is returned

    def to_rna(self): # Function converts DNA sequence object to RNA sequence
        rna_seq = '' # Initializes RNA sequence
        for base in self.sequence: # Loop goes through each nucleotide in DNA sequence object
            if base in 'ATGC':
                rna_seq += base.replace('T', 'U') # Replaces T with U
            else:
                raise ValueError('Invalid DNA sequence') # Error raised if no DNA sequence object
        return rna_seq # Returns RNA sequence

    def to_protein(self): # Function converts DNA sequence object to protein sequence
        codon_table = { #... Initializes dictionary of codons and corresponding amino acids
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '', 'UAG': '',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'UGU': 'C', 'UGC': 'C', 'UGA': '', 'UGG': 'W',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        rna_seq = self.to_rna() # Convert DNA sequence to RNA sequence
        
        protein_seq = '' #Initialize protein sequence
        for i in range(0, len(rna_seq), 3): # For loop iterates through RNA sequence in increments of 3 to generate codon sequence
            codon = rna_seq[i:i + 3] # Codon sequence is 3 nucleotides of RNA
            if len(codon) == 3:
                amino_acid = codon_table.get(codon, 'X') # Uses codon to get amino acid in codon dictionary
                protein_seq += amino_acid # Adds amino acid to amino acid sequence
        return protein_seq

def predict_protein(prot_seq): # Function sends protein sequence to SwissModel to get 3D predictive model
    token = "7fe1121e540e47550157b672eb9d17c9cfcfbe51" # API token to access SwissModel
       
    target_sequences = [prot_seq] # Define target protein sequence to be sent to SwissModel
    print("Protein sequence sent to Swiss-Model:",prot_seq)
    project_title = "Predictive 3D Protein Model" # Define project title
    
    response = requests.post( # Send request to SwissModel to start 3D protein model creation 
        "https://swissmodel.expasy.org/automodel",
        headers={ "Authorization": f"Token {token}" },
        json={ "target_sequences": target_sequences, "project_title": project_title })
    
    project_id = response.json()["project_id"] # Obtain project_id from response 
        
    import time # Loop until project complete
    while True:
        time.sleep(10)
        response = requests.get(  # Update status from server
            f"https://swissmodel.expasy.org/project/{ project_id }/models/summary/", 
            headers={ "Authorization": f"Token {token}" })
        status = response.json()["status"] # Update status
        print('Job status is now', status)
        if status in ["COMPLETED", "FAILED"]:
            break
    
    response_object = response.json() # Check if job completed
    if response_object['status'] == 'COMPLETED':
        for model in response_object['models']:
            coordinates_url = model['coordinates_url']
            
    pdb_gz_file = os.path.basename(coordinates_url) # Download PDB file
    with requests.get(coordinates_url, stream=True) as r:
        r.raise_for_status()
        with open(pdb_gz_file, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
            
    pdb_file = pdb_gz_file.replace('.gz', '') # Extract PDB file
    with gzip.open(pdb_gz_file, 'rb') as f_in: 
        with open(pdb_file, 'wb') as f_out:
            f_out.write(f_in.read())
    
    subprocess.run(['pymol', pdb_file]) # Open PDB file in PyMOL

    pymol.cmd.sync()
    pymol.cmd.quit()

def main():
    dna_seq = "" # Initialize DNA sequence variable
    while True:
        print('\nMenu Options:')
        print('1. Input a DNA sequence')
        print('2. Convert DNA to RNA sequence')
        print('3. Convert DNA to protein sequence')
        print('4. Generate a 3D protein model')
        print('5. Quit\n')
        
        choice = input('Enter your choice: ')  # Get user menu option input
        
        if choice == '1':
            dna_input = input('\nEnter the DNA sequence: ') # User inputs DNA sequence
            if set(dna_input.upper()).issubset(set("ATGC")): # Verifies DNA sequence
                dna_seq = dna_input.upper() # DNA sequence converted to uppercase
                dna = DNA(dna_seq) # DNA sequence object created
                count = dna.count_nucleotides() # Counts each nucleotide in DNA sequence
                print("\nNucleotide count: A={}, T={}, G={}, C={}".format(count[0], count[1], count[2], count[3]))
            else: # If DNA sequence is invalid, return to menu
                print('\nInvalid DNA sequence. Please try again.')

        elif choice == '2':
            if dna_seq: # Checks if DNA sequence inputed by user in choice 1
                dna = DNA(dna_seq) # DNA sequence object created
                rna_seq = dna.to_rna() # DNA converted to RNA
                print('\nRNA sequence:', rna_seq) 
            else: # If user did not input DNA sequence in choice 1, return to menu
                print("\nPlease input a DNA sequence first.")
                continue 

        elif choice == '3':
            if dna_seq: # Checks if DNA sequence inputed in choice 1
                dna = DNA(dna_seq) # DNA sequence object created
                prot_seq = dna.to_protein() # DNA converted to amino acid sequence
                print('\nProtein sequence:', prot_seq)
            else: # If user did not input DNA sequence in choice 1, return to menu
                print("\nPlease input a DNA sequence first.") 
                continue  # return to menu

        elif choice == '4':
            while True:
                print('\n1. Generate a 3D protein model from inputed DNA sequence')
                print('2. Input a protein sequence')
                print('3. Return to menu\n')
                              
                sub_choice = input('Enter your choice: ')
                print(' ')
                if sub_choice == '1':
                    if dna_seq: # Checks if DNA sequence inputed in choice 1
                        dna = DNA(dna_seq) # DNA sequence object created
                        prot_seq = dna.to_protein() # DNA converted to amino acid sequence                    
                        
                        try: # Generate a 3D model from protein sequence
                            predict_protein(prot_seq)
                            print('Protein model generated successfully.')
                        except Exception as e: # Throws exception if Swiss-Model rejects protein sequence
                            print('Error generating protein model:', e)
                        pymol.cmd.delete('*')
                    else: # Return to menu if no protein sequence
                        print('Error. Please input a protein sequence.')
                        
                        
                elif sub_choice == '2':
                    prot_seq2 = input('Input a protein sequence: ')
                    print ('')
                    try: # Generate a 3D model from inputed protein sequence
                        predict_protein(prot_seq2)
                        print('Protein model generated successfully.')
                    except Exception as e: # Throws exception if Swiss-Model rejects protein sequence
                        print('Error generating protein model: ', e)
                    
                    break                                      
                elif sub_choice == '3':
                    
                    break
                else:
                    print('Invalid choice. Please try again.')
                    continue

        elif choice == '5': # Quits program 
            break
        else: 
            print('Invalid choice. Please try again.')
main()
