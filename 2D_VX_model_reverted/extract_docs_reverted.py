import os
import glob

def extract_matlab_doc(filepath):
    doc_lines = []
    function_def = None
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('function'):
                function_def = stripped
            elif stripped.startswith('%') or stripped.startswith('%%'):
                doc_lines.append(stripped)
            elif stripped and not stripped.startswith('function') and not stripped.startswith('%'):
                if doc_lines:
                    break
    
    return function_def, '\n'.join(doc_lines)

def main():
    m_files = glob.glob('**/*.m', recursive=True)
    with open('matlab_docs_summary.txt', 'w', encoding='utf-8') as out_f:
        for fpath in sorted(m_files):
            func_def, docs = extract_matlab_doc(fpath)
            out_f.write(f"=== {fpath} ===\n")
            if func_def:
                out_f.write(f"Function: {func_def}\n")
            out_f.write(f"Docs:\n{docs}\n\n")

if __name__ == '__main__':
    main()
