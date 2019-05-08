# an_min_hist

### Script

Relying on a MDAnalysis backend, the script computes the histogram of the distances between the oxygen of water molecules (default name OW, customisable) and a group defined in any MDAnalysis compatible way (given by the selection option). The histogram is computed at each frame of the trajectory and the results averaged. Options allow to control the number of bins of the histograms and the maximum distance to analyse (to avoid spending time in the bulk water regions).

```
usage: an_mindist_hist.py [-h] --pdb PDB --trj TRJ --selection SEL
                          [--rmax RMAX] [--oxwat OXWAT] [--out OUT]
                          [--nbins OUT]

  -h, --help 		show this help message and exit
  --pdb PDB		Topology file (any accepted by MDAnalysis)
  --trj TRJ			Trajectory file (any accepted by MDAnalysis)
  --selection SEL	MDAnalysis style selection for the group to analyse
  --rmax RMAX		Maximum distance up to which analyse
  --oxwat OXWAT	Atom name for water oxygen. Default "OW"
  --out OUT		Name for txt output file with histogram. Default "average_water_histo.txt"
  --nbins OUT		Nr of bins in the histogram. Default 1200
```
  
  Requirements: numpy, MDAnalysis