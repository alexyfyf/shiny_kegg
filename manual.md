## Input file format
* Please input file in csv format
* Can be either one column csv, indicating gene symbols (eg. Lmo2)
* Or can be two column csv file, and the first column is gene symbol, second column is logFC (or any other measure for your gene list, used for coloring the network)

## Species
* Only human and mouse curently

## Supported plot format
* png
* pdf, vector graph
* svg, vector graph, recommended for Mac Adobe Illustrator

## Why no figures?
* The reason is likely to be two
* First check if your species is selected correctly
* Then you can toggle the `pvalue` and `qvalue` cutoff to 1 to check the enrichment of pathways in the `Table` panel
* Figure will only be generated after you choose proper parameters for analysis and click `Run!`

## What packages are used?
* See `SessinInfo` page for details

## Contact
* Any questions or suggestions please contact feng.yan@monash.edu