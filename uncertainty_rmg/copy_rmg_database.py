# Sevy's bastardization of https://github.com/python/cpython/blob/main/Lib/shutil.py

import os
from shutil import copy2, ignore_patterns, rmtree
import fnmatch
import time
import sys

def copy_rmg_database(
    RMG_db_folder, 
    output_path,
    perturb_dict,
    N,
    ):

    def hardcopy_patterns(*patterns):
        def _hardcopy_patterns(path, names):
            matched_names = []
            for pattern in patterns:
                matched_names.extend(fnmatch.filter(names, pattern))
            return set(matched_names)
        return _hardcopy_patterns

    # start by copying a directory but ignoring the .git
    def copytree_sym(src, dst, ignore=None, hardcopy=None, symlinks=False):
        """
        Copies a tree mostly symbolically. It will make a physical copy of the
        files specified, and not copy any of the files on the ignore list.
        """
        names = os.listdir(src)
        os.makedirs(dst)

        if ignore is not None:
            ignored_names = ignore(os.fspath(src), names)
        else:
            ignored_names = set()

        if hardcopy is not None:
            hardcopy_names = hardcopy(os.fspath(src), names)
        else:
            hardcopy_names = set()

        errors = []
        for name in names:
            if name in ignored_names:
                continue

            srcname = os.path.join(src, name)
            dstname = os.path.join(dst, name)
            try:
                if symlinks and os.path.islink(srcname):
                    linkto = os.readlink(srcname)
                    os.symlink(linkto, dstname)
                elif os.path.isdir(srcname):
                    copytree_sym(
                        srcname, 
                        dstname, 
                        symlinks=symlinks, 
                        ignore=ignore, 
                        hardcopy=hardcopy
                        )
                else:
                    # if srcname in hardcopy:
                    # the training reactions.py files are getting copied by mistake. 
                    # showing up because they have the same name as the reaction
                    # library files
                    if any perturb_dict[""]
                    if name in hardcopy_names and "/training" not in name:
                        copy2(srcname, dstname)
                    else:
                        os.symlink(srcname, dstname)
            except OSError as why:
                errors.append((srcname, dstname, str(why)))

            # catch the Error from the recursive copytree so that we can
            # continue with other files
            except Exception as err:
                errors.extend(err.args[0])

    database_src = os.path.abspath(RMG_db_folder)
    if not os.path.exists(database_src):
        raise OSError(f'Could not find source database {database_src}')

    start_time = time.time()
    database_dest = os.path.abspath(
        output_path)+ "/db_" + str(N).zfill(4)
    if os.path.exists(database_dest):
        print(f"removing {database_dest} and replacing with new one")
        rmtree(database_dest)
        # continue
        # raise OSError(f'Destination already exists: {database_dest}')

    copytree_sym(
        database_src,
        database_dest,
        symlinks=True,
        ignore=ignore_patterns(
            '.conda',
            '.git',
            '.github',
            '.gitignore',
            '.travis.yml',
            '.vscode',
            'rules_[0-9][0-9][0-9][0-9].py',
            'reactions_[0-9][0-9][0-9][0-9].py',
            'surfaceThermoPt111_[0-9][0-9][0-9][0-9].py',
            'adsorptionPt111_[0-9][0-9][0-9][0-9].py',
        ),
        # TODO see if you can delete this hardcopy option
        # only hard copy one of these 
        hardcopy=hardcopy_patterns(
            'rules.py',
            'reactions.py',
            'surfaceThermoPt111.py',
            'adsorptionPt111.py',
        )
    )
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"Copied database in {elapsed_time} seconds")

if __name__ == "__main__":

    RMG_db_folder = sys.argv[1] 
    output_path = sys.argv[2]
    N = sys.argv[3]

    copy_rmg_database(
    RMG_db_folder, 
    output_path,
    N,
    )


