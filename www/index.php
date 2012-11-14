
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->


<p> <strong> Extension of BLR to deal with censored and binary traits </strong></p>

<p>
The BLR (Bayesian Linear Regression, <a href="http://cran.r-project.org/web/packages/BLR/index.html">http://cran.r-project.org/web/packages/BLR/index.html</a>)
package of R (<a href="http://cran.r-project.org">http://cran.r-project.org</a>) implements several 
types of Bayesian regression models,  including fixed effects, Bayesian Lasso (BL, Park and Casella 2008) 
and Bayesian Ridge Regression. BLR can only handle continuous outcomes. We have produced a modified (beta) 
version of BLR  (BGLR=Bayesian Generalized Linear Regression) that extends BLR by allowing regressions for 
binary and censored outcomes. Most of the inputs, processes and outputs are as in BLR. Two supporting 
documents (<a href="BLR.pdf">BLR.pdf</a> and <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3091623/">http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3091623/</a>) provide additional 
information about that. Here we focus on describing changes in inputs, internal process and outputs introduced 
to handle binary and censored outcomes. Users that are not familiar with BLR are strongly encouraged to 
first read BLR.pdf. Future developments will be released first in 
the R-forge webpage <a href="https://r-forge.r-project.org/projects/bglr/">https://r-forge.r-project.org/projects/bglr/</a> and 
subsequently as R-packages
</p>

<p>The software and supporting documents can be downaload here <a HREF="BGLR-beta.zip">Download</a></p>

<p> <strong> Work in progress: extension of BLR to include another shrinkage methods</strong> </p>

We are currently working in an extension of BLR to include another parametric and non parametric models, for example
BayesA, BayesB, BayesCpi, Reproducing Kernel Hilbert Spaces, etc. A snapshoot of the current development can be found 
in <a href="https://r-forge.r-project.org/R/?group_id=1525">https://r-forge.r-project.org/R/?group_id=1525</a>.

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
