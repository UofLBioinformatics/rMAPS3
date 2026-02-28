param(
    [string]$Genomes = "all",
    [string]$Manifest = "scripts/genomes.manifest.tsv",
    [string]$FastaRoot = "",
    [switch]$Force,
    [switch]$KeepArchive,
    [switch]$AllowUnverified,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

if (-not $FastaRoot) {
    if ($env:RMAPS_FASTA_ROOT) {
        $FastaRoot = $env:RMAPS_FASTA_ROOT
    } else {
        $FastaRoot = "genomedata"
    }
}

function Parse-Manifest([string]$path) {
    if (-not (Test-Path $path)) {
        throw "Manifest not found: $path"
    }
    $rows = @()
    foreach ($line in Get-Content $path) {
        $trim = $line.Trim()
        if (-not $trim -or $trim.StartsWith("#")) { continue }
        $parts = $line -split "`t"
        if ($parts.Count -lt 2) { continue }
        $rows += [pscustomobject]@{
            genome = $parts[0].Trim()
            url = $parts[1].Trim()
            sha256 = if ($parts.Count -ge 3) { $parts[2].Trim().ToLower() } else { "" }
        }
    }
    return $rows
}

function Expand-GzipFile([string]$inPath, [string]$outPath) {
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $inStream = [System.IO.File]::OpenRead($inPath)
    try {
        $gzip = New-Object System.IO.Compression.GZipStream($inStream, [System.IO.Compression.CompressionMode]::Decompress)
        try {
            $outStream = [System.IO.File]::Create($outPath)
            try {
                $gzip.CopyTo($outStream)
            } finally {
                $outStream.Dispose()
            }
        } finally {
            $gzip.Dispose()
        }
    } finally {
        $inStream.Dispose()
    }
}

$rows = Parse-Manifest $Manifest
if ($rows.Count -eq 0) {
    throw "No entries found in manifest: $Manifest"
}

$requested = @()
if ($Genomes -ne "all") {
    $requested = $Genomes.Split(",") | ForEach-Object { $_.Trim() } | Where-Object { $_ }
}

if (-not (Test-Path $FastaRoot)) {
    New-Item -ItemType Directory -Force -Path $FastaRoot | Out-Null
}

foreach ($row in $rows) {
    if ($requested.Count -gt 0 -and ($requested -notcontains $row.genome)) {
        continue
    }

    $genome = $row.genome
    $url = $row.url
    $sha = $row.sha256

    if (-not $url) {
        throw "Missing URL for genome '$genome' in manifest."
    }
    if (-not $AllowUnverified -and -not $sha) {
        throw "Missing sha256 for '$genome'. Add checksum to manifest or pass -AllowUnverified."
    }

    $targetDir = Join-Path $FastaRoot $genome
    New-Item -ItemType Directory -Force -Path $targetDir | Out-Null

    $leaf = [System.IO.Path]::GetFileName($url)
    $archivePath = Join-Path $targetDir $leaf
    $finalFa = Join-Path $targetDir ($genome + ".fa")

    Write-Host "==> $genome"
    Write-Host "    URL: $url"
    Write-Host "    Target: $finalFa"

    if ($DryRun) {
        continue
    }

    if ((-not $Force) -and (Test-Path $finalFa)) {
        Write-Host "    Skip: FASTA exists"
    } else {
        if ($Force -or -not (Test-Path $archivePath)) {
            Write-Host "    Downloading..."
            Invoke-WebRequest -Uri $url -OutFile $archivePath
        } else {
            Write-Host "    Reusing existing archive"
        }

        if ($sha) {
            $got = (Get-FileHash -Algorithm SHA256 $archivePath).Hash.ToLower()
            if ($got -ne $sha) {
                throw "SHA256 mismatch for '$genome'. Expected $sha, got $got"
            }
            Write-Host "    Checksum OK"
        } else {
            Write-Warning "No checksum for '$genome' (allowed by -AllowUnverified)"
        }

        if ($archivePath.ToLower().EndsWith(".gz")) {
            Write-Host "    Decompressing..."
            Expand-GzipFile -inPath $archivePath -outPath $finalFa
        } else {
            Copy-Item -Force $archivePath $finalFa
        }
    }

    $faiPath = $finalFa + ".fai"
    if (-not (Test-Path $faiPath)) {
        $pythonCmd = Get-Command python -ErrorAction SilentlyContinue
        if ($pythonCmd) {
            Write-Host "    Indexing with pyfaidx..."
            python -c "from pyfaidx import Fasta; Fasta(r'$finalFa'); print('indexed')"
        } else {
            Write-Warning "python not found in PATH; skipping FASTA index for '$genome'"
        }
    }

    if ((-not $KeepArchive) -and (Test-Path $archivePath) -and ($archivePath -ne $finalFa)) {
        Remove-Item -Force $archivePath
    }
}

Write-Host "Done."
