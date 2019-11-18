<?PHP

/*

This program is used to plot graphs per distance. The graphs can be one per zone or five zones per graph.

*/

	require_once( dirname(__FILE__) . '/Filesystem.class.php' );

	$sSrcPath = '';
	$bGraphOnly = FALSE;
	$bSplitFilesOnly = FALSE;
	$bAllSteps = FALSE;
	$aGraphFiles = array();
	$bSplitGraphs = TRUE;

	checkRuntimeParams();

	echo "\n\nSource Dir: {$sSrcPath}\n";

	if($bAllSteps) { echo "\nMode: Split files & generate graphs"; }
	if($bGraphOnly) { echo "\nMode: Generate graphs only"; }
	if($bSplitFilesOnly) { echo "\nMode: Split files only"; }
        if($bSplitGraphs) { echo "\nMode: Splitting graphs. Multiple graphs files per distance, max five files per graph"; }
        if(!$bSplitGraphs) { echo "\nMode: Not splitting graphs. One graph file per distance"; }

	echo "\n\n\n";

	$oFS = new Filesystem();

	if($bAllSteps || $bSplitFilesOnly) {
    	
		$aFiles = $oFS->readDir($sSrcPath, '#\.dat$#');

		echo "\nSPLITTING FILES\n\n";
    
		foreach($aFiles as $sFilename) {
    
			echo "\n\tProcessing file: {$sFilename}";

			splitFile($sFilename);
    		}
	}

	if($bAllSteps || $bGraphOnly) {

		echo "\n\n\nGENERATING GRAPHS\n\n";

        	generateGraphs();
	}

    	echo "\n\nProgram ended.\n\n";



/************************************************/


function checkRuntimeParams() {

    global $argc, $argv, $sSrcPath, $sDestPath, $bGraphOnly, $bSplitFilesOnly, $bAllSteps, $bSplitGraphs;

    $bStatus = TRUE;
    
    if($argc < 2) {
	
	echo "\n\nInvalid command syntax.\n\n\nUsage: php split_nufnu_zone.php <srcfolder> [-a/A/g/G/s]\n\n-a : Split files & generate graphs. Multiple graphs files per distance, max five files per graph.\n-A : Split files & generate graphs. All graphs go into same file.\n-g : Generate graphs only. Multiple graphs files per distance, max five files per graph.\n-G : Generate graphs only. All graphs go into same file.\n-s : Default. Split files and graphs.\n\n\n";
	$bStatus = FALSE;
    }

    //echo "\n\nArg count = {$argc}\n\n";
    //print_r($argv);

    for ($i = 2; $i < $argc; $i++) {

    	echo "\n\nArgv of {$i} = {$argv[$i]}\n\n";
    
	switch($argv[$i]) {

		case '-a':
			$bAllSteps = TRUE;
			break;

		case '-A':
		        $bAllSteps = TRUE;
			$bSplitGraphs = FALSE;
			break;

		case '-g':
			$bGraphOnly = TRUE;
			break;

		case '-G':
		        $bGraphOnly = TRUE;
			$bSplitGraphs = FALSE;
                        break;

		case '-s':
		default:
     		   	$bSplitFilesOnly = TRUE;
			$bSplitGraphs = TRUE;
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
This function will read data from all the files and will put the values for each zone into one file so that the data for
each zone is in one file.
*/

function splitFile($sFile) {

	 global $sSrcPath;

	 $aLines = file("{$sSrcPath}/{$sFile}");

	 $sFilenamePrefix = '';

	 if (!preg_match('/^([^\n]+)\.dat$/', $sFile, $aMatch)) {
	 
		echo "\n\nUnable to read filename prefix.\n\n";
		exit;
	 
	 } else { 

	 	 $sFilenamePrefix = $aMatch[1];
	 }

	 

	 //print_r($aLines);

	 $iZone = FALSE;

	 foreach( $aLines as $sLine ) {

	 	if (preg_match('/^zn = ([0-9]+)\n/', $sLine, $aMatch)) {
		   
		   $iZone = $aMatch[1];

		   //echo "\nZONE = {$iZone}";

		} else if ($iZone != FALSE) {

		   file_put_contents("{$sSrcPath}/output/{$sFilenamePrefix}_zone_{$iZone}.dat", $sLine, FILE_APPEND);
		} 
	 }
}


function generateGraphs() {

	 global $sSrcPath, $aGraphFiles, $bSplitGraphs;

	 $sCurrentDir = getcwd();

         $oFS = new Filesystem();

	 $aFiles = $oFS->readDir($sSrcPath . '/output', '#\.dat$#');

         foreach($aFiles as $sFilename) {

	 	//echo "\nProcessing file: {$sFilename}";

		if( preg_match('/(nuFnucomp_fs[0-9]+)_([^_]+)_zone_([0-9]+).dat/', $sFilename, $aMatch) ) {
		    
		    //print_r($aMatch);

		    $aGraphFiles[$aMatch[2]][$aMatch[3]] = $sFilename;
		}
	}

	$sPrefix = $aMatch[1];

	//print_r($aGraphFiles);

	chdir("{$sSrcPath}/output");

	foreach ($aGraphFiles as $sZone => $aFiles) {
		
		$sScript = '';

		$iCount = 0;
		$iGraphFileCount = 1;

		//print_r($aFiles);

		ksort($aFiles);

		//print_r($aFiles);

		foreach($aFiles as $iZone => $sFile) {
			
			$iGraphCount = $iCount % 5;
			
			if ($bSplitGraphs && $iGraphCount == 4) {

				$sScript .= " -nxy {$sFile} ";

				$sCommand = "gracebat {$sScript} -log xy -autoscale xy -saveall {$sPrefix}_{$sZone}_graph_{$iGraphFileCount}.agr";

				$iGraphFileCount++;

				echo "\n\n{$sCommand}\n\n";

				echo shell_exec($sCommand);

		 	        $sScript = '';
			
			} else {

                                $sScript .= " -nxy {$sFile} ";
			}
			
			$iCount++;
		}

		if(!empty($sScript)) {

			$sCommand = "gracebat {$sScript} -log xy -autoscale xy -saveall {$sPrefix}_{$sZone}_graph_{$iGraphFileCount}.agr";

                        echo "\n\n{$sCommand}\n\n";

			echo shell_exec($sCommand);
		}

		//gracebat -nxy nuFnucomp_fs29_zht6.006692e+17_zone_96.dat -log xy -autoscale xy -savefile grace.agr
	}

	


	chdir($sCurrentDir);
}
