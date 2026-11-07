"""Recompute formulas or render saved results for Paper06."""
import argparse

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    comp=sub.add_parser('recompute',help='Derive exact results into a new directory')
    comp.add_argument('--output-dir',required=True)
    fig=sub.add_parser('render',help='Render supplied data into a new directory')
    fig.add_argument('--data-dir',required=True);fig.add_argument('--output-dir',required=True)
    args=parser.parse_args()
    if args.command=='recompute':
        from paper06_results.results import recompute
        meta=recompute(args.output_dir)
        print(f"Completed {meta['identities']} derivation checks (including numerical controls) and saved {len(meta['files'])} data files.")
    else:
        from paper06_results.render import render
        render(args.data_dir,args.output_dir);print('Rendered six figures from saved data.')

if __name__=='__main__':main()
