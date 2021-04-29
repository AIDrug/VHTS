#!/usr/bin/env python
import sys
import argparse
from vhts.pydock import cal_box_size


class LoadFromConfig(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with values as f:
            parser.parse_args(f.read().split(), namespace)


class ExtendAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)


def main():

    parser = argparse.ArgumentParser(
        description='cal docking box from ligand files')
    parser.register('action', 'extend', ExtendAction)
    parser.add_argument('-i', '--autobox_file_list', type=str, action='extend',
                        nargs='+', required=False, default=[],
                        help='--autobox_file_list a.pdb b.pdb\n' +
                             '--autobox_file_list c.pdb')
    parser.add_argument('-n', '--autobox_margin', type=float, required=False,
                        default=4.0, help='autobox margin default: 4.0')
    parser.add_argument('--arg_file', type=open, required=False, default=None,
                        action=LoadFromConfig, help='argment file')

    args = parser.parse_args()
    autobox_file_list = args.autobox_file_list
    autobox_margin = args.autobox_margin
    if len(autobox_file_list) < 1:
        print('the following arguments are required: --autobox_file_list')
        sys.exit()
    box_center, box_size = cal_box_size(
        autobox_file_list, margin=autobox_margin, use_hydrogen=False)
#    box_parameter = (box_center, box_size)

    line_out = 'center_x=%.3f\ncenter_y=%.3f\ncenter_z=%.3f\n' % (box_center)
    line_out += 'size_x=%.3f\nsize_y=%.3f\nsize_z=%.3f\n' % (box_size)

    print(line_out.strip())


if __name__ == "__main__":
    main()
