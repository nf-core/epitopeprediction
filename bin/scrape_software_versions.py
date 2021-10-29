#!/usr/bin/env python
from __future__ import print_function
import os

regexes = {
    'nf-core/epitopeprediction': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'CSVTK': ['v_csvtk.txt', r"csvtk v(\S+)"],
    'SNPsift': ['v_snpsift.txt', r"SnpSift version (\S+)"],
    'Fred2': ['v_fred2.txt', r"fred2 (\S+)"],
    'MHCFlurry': ['v_mhcflurry.txt', r"mhcflurry (\S+)"],
    'MHCnuggets': ['v_mhcnuggets.txt', r"mhcnuggets (\S+)"],
    'NetMHC': ['v_netmhc.txt', r"netmhc (\S+)"],
    'NetMHCpan': ['v_netmhcpan.txt', r"netmhcpan (\S+)"],
    'NetMHCII': ['v_netmhcii.txt', r"netmhcii (\S+)"],
    'NetMHCIIpan': ['v_netmhciipan.txt', r"netmhciipan (\S+)"],
}

results = OrderedDict()
results['nf-core/epitopeprediction'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['CSVTK'] = '<span style="color:#999999;\">N/A</span>'
results['SNPsift'] = '<span style="color:#999999;\">N/A</span>'
results['Fred2'] = '<span style="color:#999999;\">N/A</span>'
results['MHCFlurry'] = '<span style="color:#999999;\">N/A</span>'
results['MHCnuggets'] = '<span style="color:#999999;\">N/A</span>'
results['NetMHC'] = '<span style="color:#999999;\">N/A</span>'
results['NetMHCpan'] = '<span style="color:#999999;\">N/A</span>'
results['NetMHCII'] = '<span style="color:#999999;\">N/A</span>'
results['NetMHCIIpan'] = '<span style="color:#999999;\">N/A</span>'

version_files = [x for x in os.listdir(".") if x.endswith(".version.txt")]
for version_file in version_files:
    software = version_file.replace(".version.txt", "")
    if software == "pipeline":
        software = "nf-core/epitopeprediction"

    with open(version_file) as fin:
        version = fin.read().strip()
    results[software] = version

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/epitopeprediction Software Versions'
section_href: 'https://github.com/nf-core/epitopeprediction'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in sorted(results.items()):
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out as tsv file:
with open("software_versions.tsv", "w") as f:
    for k, v in sorted(results.items()):
        f.write("{}\t{}\n".format(k, v))
