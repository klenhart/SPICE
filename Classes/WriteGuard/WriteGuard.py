#!/bin/env python

#######################################################################
# Copyright (C) 2023 Christian Bluemel
#
# This file is part of Spice.
#
#  WriteGuard is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  WriteGuard is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Spice.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

import os
import time
import errno

from Classes.WriteGuard.NoHashBut import NoHashBut


class WriteGuard:

    def __init__(self, guarded_file_path: str, guarded_file_dir_path: str, tag: str = ""):
        self.tag = tag
        self.guarded_flag = False
        self.guarded_file_path = guarded_file_path
        no_hash: NoHashBut = NoHashBut(guarded_file_path)
        self.guard_tag = no_hash.represent
        self.guard_temp_file_name = os.path.join(guarded_file_dir_path, "guard" + str(self.guard_tag) + ".temp")

    def end_guard(self) -> None:
        self.guarded_flag = False
        os.remove(self.guard_temp_file_name)

    def start_guard(self) -> bool:
        start_time: float = time.time()
        time_limit: int = 3600
        sleep_time: float = 31.5
        sleep_deduction: float = 1.0
        # Run this until the guard could be started.
        while True:
            try:
                temp_guard_file = os.open(self.guard_temp_file_name, os.O_CREAT | os.O_EXCL | os.O_RDWR)
                with os.fdopen(temp_guard_file, "a") as f:
                    f.write(self.guarded_file_path)
                break
            except OSError as oe:
                if oe.errno != errno.EEXIST:
                    raise
                if (time.time() - start_time) >= time_limit:
                    print("Establishing lock timed out. Took longer than", time_limit, "seconds.")
                    return False
                print(self.tag, "Guard already in place. Passed time:", time.time() - start_time)
                time.sleep(sleep_time)
                # The longer a process has to wait, the faster the request frequency becomes.
                if sleep_time > 0.5:
                    sleep_time = sleep_time - sleep_deduction
                    sleep_deduction = sleep_deduction*2
        self.guarded_flag = True
        return True

    def __enter__(self):
        self.start_guard()

    def __exit__(self, ignore, value, traceback):
        self.end_guard()


if __name__ == "__main__":
    directory: str = "some\\path\\"
    file: str = "some\\path\\"
    with WriteGuard(file, directory):
        writeGuardAlt = WriteGuard(file, directory)
        writeGuardAlt.start_guard()
