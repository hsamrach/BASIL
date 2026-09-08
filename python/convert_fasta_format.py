#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path


def convert_fasta_file(input_file, output_file):

    try:
        with open(input_file, 'r') as infile:
            content = infile.read()
        
        # should start with '>'
        if not content.strip().startswith('>'):
            print(f"ERROR: {input_file} does not appear to be a valid FASTA file (no '>' header)")
            return False
        
        # Write to output
        with open(output_file, 'w') as outfile:
            outfile.write(content)
        
        return True
    
    except IOError as e:
        print(f"ERROR: Failed to process {input_file}: {e}")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error processing {input_file}: {e}")
        return False


def convert_directory(input_dir, output_dir):
    
    #Convert all FASTA files in a directory from non-.fasta formats to .fasta format.
    supported_extensions = {'.fas', '.fa', '.fna'}
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    
    if not input_path.is_dir():
        print(f"ERROR: Input directory does not exist: {input_dir}")
        return {'total': 0, 'converted': 0, 'skipped': 0, 'failed': 0, 'conversions': []}
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    summary = {'total': 0, 'converted': 0, 'skipped': 0, 'failed': 0, 'conversions': []}
    
    # Find all FASTA files
    for file_path in sorted(input_path.iterdir()):
        if not file_path.is_file():
            continue
        
        suffix = file_path.suffix.lower()
        stem = file_path.stem
        
        summary['total'] += 1
        
        # If already .fasta, skip
        if suffix == '.fasta':
            output_file = output_path / file_path.name
            try:
                # Still copy to maintain consistency
                with open(file_path, 'r') as infile:
                    content = infile.read()
                with open(output_file, 'w') as outfile:
                    outfile.write(content)
                summary['skipped'] += 1
                summary['conversions'].append({
                    'input': file_path.name,
                    'output': output_file.name,
                    'status': 'skipped (already .fasta)'
                })
            except Exception as e:
                print(f"ERROR: Failed to copy {file_path.name}: {e}")
                summary['failed'] += 1
                summary['conversions'].append({
                    'input': file_path.name,
                    'output': None,
                    'status': f'failed ({e})'
                })
            continue
        
        # If supported extension, convert
        if suffix in supported_extensions:
            output_file = output_path / f"{stem}.fasta"
            if convert_fasta_file(str(file_path), str(output_file)):
                summary['converted'] += 1
                summary['conversions'].append({
                    'input': file_path.name,
                    'output': output_file.name,
                    'status': 'converted'
                })
            else:
                summary['failed'] += 1
                summary['conversions'].append({
                    'input': file_path.name,
                    'output': None,
                    'status': 'failed'
                })
        else:
            # Ignore files with unsupported extensions
            summary['total'] -= 1
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Convert FASTA files from .fas, .fa, .fna to .fasta format'
    )
    parser.add_argument(
        'input_dir',
        help='Directory containing FASTA files to convert'
    )
    parser.add_argument(
        'output_dir',
        help='Directory to write converted .fasta files'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed conversion information'
    )
    
    args = parser.parse_args()
    
    print(f"Converting FASTA files from {args.input_dir} to {args.output_dir}")
    summary = convert_directory(args.input_dir, args.output_dir)
    
    # Print summary
    print("\n" + "="*60)
    print("CONVERSION SUMMARY")
    print("="*60)
    print(f"Total files processed:  {summary['total']}")
    print(f"Files converted:        {summary['converted']}")
    print(f"Files skipped:          {summary['skipped']}")
    print(f"Files failed:           {summary['failed']}")
    
    if args.verbose and summary['conversions']:
        print("\nDetailed conversion log:")
        print("-"*60)
        for conv in summary['conversions']:
            status = conv['status']
            input_name = conv['input']
            output_name = conv['output'] if conv['output'] else 'N/A'
            print(f"  {input_name:40} → {output_name:40} [{status}]")
    
    if summary['failed'] > 0:
        print(f"\nWARNING: {summary['failed']} file(s) failed to convert")
        sys.exit(1)
    
    if summary['converted'] + summary['skipped'] == 0:
        print("WARNING: No FASTA files to convert")
        sys.exit(1)
    
    print("\nConversion completed successfully!")
    sys.exit(0)


if __name__ == '__main__':
    main()
