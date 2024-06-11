<p align="justify"><b>This tutorial is designed to calculate the energy contribution of amino acid residues to the energy barrier of a given reaction. This can be achieved by performing single-point calculations on the provided reactant and transition state structures, each with the specified residue deleted. </b></p>

<p align="justify"> It requires a file with a list of residues to be deleted, a *prmtop file, reactant and TS structures in the *pdb format, a CP2K input template and a file with a VMD selection of the QM region. The following packages are also required: VMD, cpptraj, parmed. </p>

---

<br>
<h2> <p align="center"> <b>I - Input Preparation </b> </p></h2>

<br/>

The <a href="https://arvpinto.github.io/enzyme_ts_macrodipole_cp2k/res_qmmm_cp2k.sh" target="_blank">res_qmmm_cp2k.sh</a> script has the following usage:

```js
./res_qmmm_cp2k.sh res_list.dat topology.prmtop R.pdb TS.pdb cp2k_template.inp qm_selection.dat
```
It prepares a directory for each residue in the list where the input files for CP2K will be output. The supplied topology and structures will be processed through cpptraj to delete each of the specified residues. Avoid deleting residues from the QM/MM boundary and QM layer. Since deleting residues changes the atom numbering, the QM/MM settings have to be updated for each deletion. The <a href="https://arvpinto.github.io/enzyme_ts_macrodipole_cp2k/vmd_cp2k-qmmm.tcl" target="_blank">vmd_cp2k-qmmm.tcl</a> script is called within the latter to produce a file with the configuration of the QM layer, defined by the selection in the qm_selection.dat file. It also prepares a parmed input file that corrects the charges in the *prmtop file for the electrostatic embedding scheme. The cp2k_template.inp file should have tags (PRMTOP_TAG and STATE_TAG) that are replaced with the corresponding filenames. 

<br/>

<h2> <p align="center"> <b>II - Output Processing</b> </p></h2>

<br/>

After running the single-point calculations, the following command allows to extract the absolute energies and calculate the R->TS energy barrier for each residue deletion:

```js
paste <(for i in RES_*; do echo "$i" | sed 's/RES_//g'; done) <(for i in RES_*; do echo $(grep "Total FORCE" "$i"/res_qmmm_TS.out | tail -n -1) ; done | awk '{print $9}') <(for i in RES_*; do echo $(grep "Total FORCE" "$i"/res_qmmm_R.out | tail -n -1) ; done | awk '{print $9}') | awk '{print $1,($2-$3)*627.509-14.8}' | sort -n -k1,1 > energy_differences.dat
```

<br/>

The energy barriers can be ploted with the <a href="https://arvpinto.github.io/enzyme_ts_macrodipole_cp2k/E_diff_bar_plot.py" target="_blank">E_diff_bar_plot.py</a> script:

```js
python E_diff_bar_plot.py energy_differences.dat
```

<br/>

<p align="justify"> should be the processed pc.pdb file, &lt;eps&gt; and &lt;min_samples&gt; define the parameters for outlier identification using the DBSCAN method, and &lt;n_components&gt; defines the number of clusters in the Gaussian Mixture Models clustering. The script produces a 3D plot of the PCA vectors, where the outliers are represented as black markers, the frames closest to the highest density points as white markers, and each cluster displays a different color. Additionally, the density distribution curves of each cluster are plotted against each PCA vector, with markers representing the identified frames.
Initially try different &lt;eps&gt; and &lt;min_samples&gt; values to see which and how many frames are being identified as outliers.
Once you have an adequate number of outliers, try different &lt;n_components&gt; values to identify which number of clusters is more suitable.
Also take a look at the kernel density plots to see if the density distributions have a regular shape, and the identified frames lie close to highest density points. </p>
<br/>

<div align="center">
    <img src="bar_plot.png">
</div>
<br/>

<div align="center">
    <img src="marker_plot.png">
</div>
<br/>

```js
Number of DBSCAN outliers: 29
Total number of clusters (GMM): 4
Cluster 0: 595 frames
Top 5 closest frames for Cluster 0: [ 578  721  681  647 1544]
Cluster 1: 1198 frames
Top 5 closest frames for Cluster 1: [1232 1380 1293 1919 1708]
Cluster 2: 463 frames
Top 5 closest frames for Cluster 2: [114  69  68  64  67]
Cluster 3: 215 frames
Top 5 closest frames for Cluster 3: [2015 2076 2050 2052 2054]
```

<br>

<br>
<br/>

A clusters.csv file is outputed with the cluster numbers that each frame corresponds to (outliers belong in the -1 cluster).
A frames.dat is ouputed with the top 5 frames that are closest to the highest density point of each cluster.

<br>
<h2> <p align="center"> <b>III - Frame extraction</b> </p></h2>

<br/>

Use the <a href="https://arvpinto.github.io/3D_clustering_PCA/extract_highdens.py" target="_blank">extract_highdens.py</a> script to extract the identified frames from the trajectory.
The <a href="https://arvpinto.github.io/3D_clustering_PCA/extract_highdens.py" target="_blank">extract_highdens.py</a> script usage follows:

```js
python extract_highdens.py <xtc_file> <gro_file> <cluster_indices_file> <output_prefix>
```






