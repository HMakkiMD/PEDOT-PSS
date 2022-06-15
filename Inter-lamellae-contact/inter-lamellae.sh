#!/usr/bin/expect -f

set timeout -1
send -- "export GMX_MAXBACKUP=-1"

expect "$"
spawn mkdir GRO
expect "$"
spawn mkdir STACK

set counter 0
set trajframes 10
set frame 25012

while { $counter < $trajframes } {
expect "$"
spawn gmx trjconv -f traj.xtc -s traj.tpr -o GRO/$counter.gro -b $frame -e $frame -pbc mol

expect "Select a group:"
send -- "0\r"

expect "$ "
spawn python3 inter-lamellae.py
expect "Enter EDOT heavy atoms index file name: "
send -- "pedot-heavy-all.ndx\r"
expect "Enter EDOT sp2 carbon atoms index file name: "
send -- "sp2c-all.ndx\r"
expect "Enter a gro file name: "
send -- "GRO/$counter.gro\r"


expect "$ "
spawn mv cluster.ndx STACK/cluster$counter.ndx
expect "$ "
spawn mv isolatedchains.ndx STACK/isolatedchains$counter.ndx
expect "$ "
spawn mv network.txt STACK/network$counter.ndx
expect "$ "
spawn mv pistack.ndx STACK/pistack$counter.ndx
expect "$ "
spawn mv pistack.txt STACK/pistack$counter.txt
expect "$ "
spawn mv information.txt STACK/information$counter.txt
expect "$ "
spawn mv NUMBEROFCONTACT.txt STACK/NUMBEROFCONTACT$counter.txt
expect "$ "
spawn mv LARGESTCLUSTER.txt STACK/LARGESTCLUSTER$counter.txt
expect "$ "
spawn mv CONNECTEDNETWORK.txt STACK/CONNECTEDNETWORK$counter.txt
expect "$ "
spawn rm cluster.ndx isolatedchains.ndx network.txt pistack.ndx pistack.txt information.txt NUMBEROFCONTACT.txt LARGESTCLUSTER.txt CONNECTEDNETWORK.txt

set counter [ expr $counter+1 ]
set frame [ expr $frame+104 ]
}

expect eof