<?PHP

/*

This programs will generate graphs from nuFnucomp*.dat files. 
You can generate one graph per file or generate one graph for all the files in the source directory.
The graph will plot all five data columns in the .dat file. Select columns can be ignored when plotting by setting the right options.

*/

	require_once( dirname(__FILE__) . '/Filesystem.class.php' );
	
	$sScriptFolder = dirname(__FILE__);
	$sSrcPath = '';
	$bMultiFileGraph = FALSE;
	$bIgnoreColumns = FALSE;

	checkRuntimeParams();

	echo "\n\nSource Dir: {$sSrcPath}\n";

        if(!$bMultiFileGraph) { echo "\nMode: Generate one graph per dat file."; }
	if($bMultiFileGraph) { echo "\nMode: Generate one graph for all dat files."; }
	if($bIgnoreColumns) { echo "\nMode: Ignoring specified columns in the dat file."; }
        if(!$bIgnoreColumns) { echo "\nMode: Plotting all five columns of data in the dat files."; }

	echo "\n\n\n";

	if($bMultiFileGraph == FALSE) {

		$oFS = new Filesystem();

		$aFiles = $oFS->readDir($sSrcPath, '#\.dat$#');

			echo "\nPLOTTING FILES\n\n";
    
		foreach($aFiles as $sFilename) {
    
			plotFile($sFilename);
		}

	} else {

	       echo "\nPLOTTING FOLDER\n";

	       plotFolder($sSrcPath);
	}

    	echo "\n\nProgram ended.\n\n";



/************************************************/


function checkRuntimeParams() {

    global $argc, $argv, $sSrcPath, $bMultiFileGraph, $bIgnoreColumns;

    $bStatus = TRUE;
    
    if($argc < 2) {
	
	echo "\n\nInvalid command syntax.\n\nUsage: php generate_nufnucomp_graphs [-m | -i[1-5] ] <source folder>\n";
	echo "\n-m : Generate one graph for all data files data.";
	echo "\n-i[1-5] : Ignore column specified the number in the dat files. Example: -i3 will ignore third column data. -i2 -i5 will ignore second and fifth columns in the data.";
	$bStatus = FALSE;
    }

    //echo "\n\nArg count = {$argc}\n\n";
    //print_r($argv);

    for ($i = 2; $i < $argc; $i++) {

    	//echo "\n\nArgv of {$i} = {$argv[$i]}\n\n";
    
	switch($argv[$i]) {

		case '-m':
			$bMultiFileGraph = TRUE;
			break;

		case '-i':
			$bIgnoreColumns = TRUE;
			break;

		default:
                        $bMultiFileGraph = FALSE;
	                $bIgnoreColumns = FALSE;
		     	break;
	}
    }
    
    $sSrcPath = realpath($argv[1]);
    
    if( !is_dir($sSrcPath) ) {
	
	echo "\nError: Source folder error. Please check \"{$sSrcPath}\" is valid and exists.";
	$bStatus = FALSE;
    }

    
    if( is_dir($sSrcPath) && !is_dir("{$sSrcPath}/output/") && mkdir("{$sSrcPath}/output/") == FALSE) {
    
	$bStatus = FALSE;
    }
    
    if($bStatus !== TRUE) {
	
	echo "\n\nPlease fix the errors above and try again.\n\n\n";
	exit;
    }
}

/*
This function will plot each file into a graph
*/

function plotFile($sFile) {

	 global $sSrcPath, $sScriptFolder;

 	echo "\nProcessing file: {$sFile}\n";

	// nuFnucomp_fs9_zht1.099568e+17.dat

	if( preg_match('/nuFnucomp_fs([0-9]+)_zht([^_]+).dat/', $sFile, $aMatch) ) {

	    //print_r($aMatch);

	    $iRunId = $aMatch[1];
	    $iDistance = $aMatch[2];
	
	} else {

	   echo "\nSkipping {$sFile}. Unable to parse file name based on expected pattern. Please check the file name.\n";
	   return;
	}

	$sCommand = "gracebat -nxy {$sFile} -param {$sScriptFolder}/graph_params.par -world 1e9 1e10 1e27 1e14 -noprint -saveall output/nuFnucomp_fs{$iRunId}_zht{$iDistance}_graph.agr";

	
	chdir("{$sSrcPath}");
	$sOutput = shell_exec($sCommand);

	//echo "\n\n===========\n{$sOutput}===========\n\n";

}


function plotFolder($sFolder) {

	 global $sSrcPath, $sScriptFolder;

 	echo "\nProcessing folder: {$sFolder}\n";

	$oFS = new Filesystem();

	$aFiles = $oFS->readDir($sFolder, '#\.dat$#');
	$aSortedFiles = array();
	foreach($aFiles as $sFile) {

		preg_match('/nuFnucomp_fs([0-9]+)_zht([^_]+).dat/', $sFile, $aMatch);
		$iKey = $aMatch[2];
		$aSortedFiles[$iKey] = $sFile;
	}
	ksort($aSortedFiles);
	$aFiles = $aSortedFiles;

	echo "\nFiles in folder:\n";
	print_r($aFiles);
	echo "\n\n";

	// nuFnucomp_fs9_zht1.099568e+17.dat

	if( preg_match('/nuFnucomp_fs([0-9]+)_zht([^_]+).dat/', $aFiles[key($aFiles)], $aMatch) ) {

	    //print_r($aMatch);

	    $iRunId = $aMatch[1];
	    $iDistance = $aMatch[2];
	
	} else {

	   echo "\nSkipping {$sFolder}. Unable to parse file name based on expected pattern. Please check the file name.\n";
	   return;
	}

	$sAllFiles = '';

	foreach($aFiles as $sFilename) {
	
		$sAllFiles .= " -nxy {$sFilename} ";
	}

	$sCommand = "gracebat {$sAllFiles} -param {$sScriptFolder}/graph_params.par -world 1e9 1e10 1e27 1e14 -noprint -saveall nuFnucomp_fs{$iRunId}_zhtALL_graph.agr";

	echo "Executing command: \n\n\t{$sCommand}\n\n";

	chdir("{$sFolder}");
	$sOutput = shell_exec($sCommand);

	//echo "\n\n===========\n{$sOutput}===========\n\n";

}

