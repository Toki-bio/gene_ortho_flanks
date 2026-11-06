#!/bin/bash
# IGV package per gene item: both haplotype windows, reads merged, exons, N run, fill
set -uo pipefail
. "$(dirname "$(readlink -f "$0")")/cfg.sh"
cd "$W/out" || exit 1
for i in $(awk -F'\t' '$2~/gene|control/{print $1}' ../items.tsv); do
	cd $i || exit 1
	rm -rf igv; mkdir igv
	cat view.fa other.fa > igv/ref.fa; samtools faidx igv/ref.fa
	samtools merge -f igv/reads.bam view.bam other.bam; samtools index igv/reads.bam
	cat view.bed other.bed > igv/ref.bed
	loc=$(awk '{printf "%s%s:1-%s", (NR>1?" ":""), $1, $2}' igv/ref.fa.fai)
	cat > igv/session.xml <<E
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="ref.fa" hasGeneTrack="false" hasSequenceTrack="true" locus="$loc" version="8">
    <Resources>
        <Resource path="reads.bam" type="bam"/>
        <Resource path="ref.bed" type="bed"/>
    </Resources>
    <Panel height="600" name="DataPanel" width="1400">
        <Track attributeKey="reads.bam Coverage" autoScale="true" id="reads.bam_coverage" name="HiFi coverage" showReference="false" visible="true"/>
        <Track attributeKey="reads.bam" id="reads.bam" name="HiFi reads ERR16751478, each window mapped alone (MAPQ>=20, aligned >=1 kb)" visible="true"/>
    </Panel>
    <Panel height="60" name="FeaturePanel" width="1400">
        <Track attributeKey="ref.bed" id="ref.bed" name="exons, N run, fill" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.8"/>
</Session>
E
	tar czf igv.tgz -C igv . ; rm -rf igv
	cd ..
done
ls -la */igv.tgz | awk '{s+=$5} END{print NR, "igv packages", s}'
