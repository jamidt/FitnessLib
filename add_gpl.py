import os

gpl = "/**\n\
 * Copyright © 2017 Jan Schmidt\n\
 * \n\
 * This program is free software: you can redistribute it and/or modify\n\
 * it under the terms of the GNU General Public License as published by\n\
 * the Free Software Foundation, either version 3 of the License, or\n\
 * (at your option) any later version.\n\
 * \n\
 * This program is distributed in the hope that it will be useful,\n\
 * but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
 * GNU General Public License for more details.\n\
 * \n\
 * You should have received a copy of the GNU General Public License\n\
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\
 */\n"

for root, dirs, files in os.walk('./'):
    if '.git' in root.split('/'):
        continue
    else:
        for filename in files:
            if filename[-4:] == ".hpp" or filename[-4:] == ".cpp":
                with open(root + "/" + filename, "r+") as f:
                    content = f.read()
                    f.seek(0, 0)
                    f.write(gpl + "\n" + content)
