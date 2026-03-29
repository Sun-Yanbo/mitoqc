"""
Unit tests for MitoQC v0.2 core modules.
Run with: pytest tests/ -v
"""
from __future__ import annotations
import pytest
from Bio.Seq import Seq

# ── Shared helpers ────────────────────────────────────────────────────────────

def _make_rec(seq="ATGC"*4000, taxon_group="Mammalia",
              organism="Homo sapiens", features=None):
    return {
        "accession": "TEST001", "organism": organism,
        "lineage": f"Eukaryota; Metazoa; Chordata; {taxon_group}",
        "taxon_group": taxon_group, "seq_len": len(seq),
        "seq": seq.upper(), "date": "01-JAN-2024",
        "platform": "Illumina", "features": features or [], "n_features": 0,
    }

def _cds(gene, start, end, strand=1):
    return {"type":"CDS","gene":gene,"start":start,"end":end,"strand":strand}

def _rrna(gene):
    return {"type":"rRNA","gene":gene,"start":0,"end":100,"strand":1}

def _trna():
    return {"type":"tRNA","gene":"trnA","start":0,"end":70,"strand":1}

def _full_features():
    pcgs = ["nad1","nad2","nad3","nad4","nad4l","nad5","nad6",
            "cox1","cox2","cox3","atp6","atp8","cytb"]
    feats  = [_cds(g, i*1000, i*1000+900) for i,g in enumerate(pcgs)]
    feats += [_rrna("12S rRNA"), _rrna("16S rRNA")]
    feats += [_trna() for _ in range(22)]
    return feats

def _clean_cds(length=1542):
    n_codons = (length - 6) // 3
    return "ATG" + "AAA" * n_codons + "TAA"


# ── M1 tests ──────────────────────────────────────────────────────────────────

class TestM1:
    from mitoqc.modules import m1_integrity as m1

    def test_pass_normal_mammalia(self):
        from mitoqc.modules import m1_integrity as m1
        rec = _make_rec(seq="ATGC"*4000, taxon_group="Mammalia")
        r = m1.run(rec)
        assert r["m1_pass"] is True

    def test_too_short(self):
        from mitoqc.modules import m1_integrity as m1
        rec = _make_rec(seq="ATGC"*1000, taxon_group="Mammalia")  # 4000 bp < 15000
        r = m1.run(rec)
        assert r["m1_pass"] is False
        assert "TOO_SHORT" in r["m1_flags"]

    def test_high_n_content(self):
        from mitoqc.modules import m1_integrity as m1
        seq = "N"*200 + "ATGC"*4000
        rec = _make_rec(seq=seq, taxon_group="Mammalia")
        r = m1.run(rec)
        assert r["m1_pass"] is False
        assert "HIGH_N" in r["m1_flags"]

    def test_gc_outlier(self):
        from mitoqc.modules import m1_integrity as m1
        seq = "GC"*8000  # 100% GC
        rec = _make_rec(seq=seq, taxon_group="Mammalia")
        r = m1.run(rec)
        assert "HIGH_GC" in r["m1_flags"]


# ── M2 tests ──────────────────────────────────────────────────────────────────

class TestM2:
    def test_pass_complete_annotation(self):
        from mitoqc.modules import m2_annotation as m2
        rec = _make_rec(features=_full_features())
        r = m2.run(rec)
        assert r["m2_pass"] is True

    def test_missing_pcg(self):
        from mitoqc.modules import m2_annotation as m2
        feats = [f for f in _full_features()
                 if not (f["type"]=="CDS" and f["gene"]=="nad1")]
        rec = _make_rec(features=feats)
        r = m2.run(rec)
        assert r["m2_pass"] is False
        assert "nad1" in r["m2_missing_pcgs"]

    def test_nematoda_atp8_absent_is_ok(self):
        from mitoqc.modules import m2_annotation as m2
        feats = [f for f in _full_features()
                 if not (f["type"]=="CDS" and f["gene"]=="atp8")]
        rec = _make_rec(features=feats, taxon_group="Nematoda",
                        organism="Caenorhabditis elegans")
        r = m2.run(rec)
        assert "atp8" not in r["m2_missing_pcgs"]

    def test_low_trna(self):
        from mitoqc.modules import m2_annotation as m2
        feats = [f for f in _full_features() if f["type"] != "tRNA"]
        feats += [_trna() for _ in range(5)]  # only 5 tRNAs
        rec = _make_rec(features=feats)
        r = m2.run(rec)
        assert r["m2_pass"] is False
        assert "LOW_TRNA" in r["m2_flags"]

    def test_missing_rrna(self):
        from mitoqc.modules import m2_annotation as m2
        feats = [f for f in _full_features() if f["type"] != "rRNA"]
        rec = _make_rec(features=feats)
        r = m2.run(rec)
        assert r["m2_pass"] is False
        assert "MISSING_RRNA" in r["m2_flags"]


# ── M3 tests ──────────────────────────────────────────────────────────────────

class TestM3:
    def test_pass_uniform_gc(self):
        from mitoqc.modules import m3_chimera as m3
        seq = "ATGC" * 4000  # uniform GC
        rec = _make_rec(seq=seq, taxon_group="Mammalia")
        r = m3.run(rec)
        assert r["m3_pass"] is True

    def test_chimera_detected(self):
        from mitoqc.modules import m3_chimera as m3
        # First half: high GC; second half: low GC
        seq = "GC"*4000 + "AT"*4000
        rec = _make_rec(seq=seq, taxon_group="Mammalia")
        r = m3.run(rec)
        assert r["m3_pass"] is False
        assert "HIGH_GC_SD" in r["m3_flags"]

    def test_insufficient_windows(self):
        from mitoqc.modules import m3_chimera as m3
        seq = "ATGC" * 200  # 800 bp — too short for 5 windows
        rec = _make_rec(seq=seq, taxon_group="Mammalia")
        r = m3.run(rec)
        assert r["m3_pass"] is True
        assert "INSUFFICIENT_WINDOWS" in r["m3_flags"]

    def test_config_window_size(self):
        from mitoqc.modules import m3_chimera as m3
        seq = "ATGC" * 4000
        rec = _make_rec(seq=seq, taxon_group="Mammalia")
        r = m3.run(rec, config={"window_size": 1000})
        assert r["m3_pass"] is True


# ── M4 tests ──────────────────────────────────────────────────────────────────

class TestM4:
    def _cox1_feat(self, start=0, length=1542):
        return _cds("cox1", start, start+length)

    def test_pass_good_record(self):
        from mitoqc.modules import m4_species as m4
        rec = _make_rec(organism="Homo sapiens",
                        features=[self._cox1_feat()])
        r = m4.run(rec)
        assert r["m4_pass"] is True

    def test_coi_alias_coi(self):
        from mitoqc.modules import m4_species as m4
        feat = _cds("COI", 0, 1542)
        rec = _make_rec(organism="Gallus gallus", features=[feat])
        r = m4.run(rec)
        assert r["m4_coi_found"] is True

    def test_coi_not_found(self):
        from mitoqc.modules import m4_species as m4
        rec = _make_rec(organism="Rattus norvegicus", features=[])
        r = m4.run(rec)
        assert r["m4_coi_found"] is False
        assert "COI_NOT_EXTRACTED" in r["m4_flags"]

    def test_ambiguous_taxon_sp(self):
        from mitoqc.modules import m4_species as m4
        rec = _make_rec(organism="Drosophila sp. 123",
                        features=[self._cox1_feat()])
        r = m4.run(rec)
        assert r["m4_pass"] is False
        assert "AMBIGUOUS_TAXON" in r["m4_flags"]

    def test_non_binomial_name(self):
        from mitoqc.modules import m4_species as m4
        rec = _make_rec(organism="Drosophila",
                        features=[self._cox1_feat()])
        r = m4.run(rec)
        assert r["m4_pass"] is False

    def test_numeric_strain(self):
        from mitoqc.modules import m4_species as m4
        rec = _make_rec(organism="Mus musculus 12345",
                        features=[self._cox1_feat()])
        r = m4.run(rec)
        assert r["m4_pass"] is False
        assert "NUMERIC_STRAIN_ID" in r["m4_flags"]


# ── M5 tests ──────────────────────────────────────────────────────────────────

class TestM5:
    def _make_cds_rec(self, cds_seq, taxon_group="Mammalia", organism="Homo sapiens"):
        seq = "A"*1000 + cds_seq + "A"*(20000 - 1000 - len(cds_seq))
        feat = _cds("cox1", 1000, 1000+len(cds_seq))
        return {
            "accession": "TEST", "organism": organism,
            "lineage": f"Eukaryota; Metazoa; Chordata; {taxon_group}",
            "taxon_group": taxon_group, "seq_len": len(seq), "seq": seq,
            "date": "01-JAN-2024", "platform": "Illumina",
            "features": [feat], "n_features": 1,
        }

    def test_pass_no_internal_stop(self):
        from mitoqc.modules import m5_pcg as m5
        cds = _clean_cds(900)
        rec = self._make_cds_rec(cds)
        r = m5.run(rec)
        assert r["m5_pass"] is True

    def test_internal_stop_detected(self):
        from mitoqc.modules import m5_pcg as m5
        cds = "ATG" + "AAA"*4 + "TAA" + "AAA"*10 + "TAA"  # internal stop
        rec = self._make_cds_rec(cds)
        r = m5.run(rec)
        assert r["m5_pass"] is False
        assert "INTERNAL_STOP" in r["m5_flags"]

    def test_no_cds_features(self):
        from mitoqc.modules import m5_pcg as m5
        rec = _make_rec(features=[])
        r = m5.run(rec)
        assert r["m5_pass"] is True
        assert r["m5_n_cds_checked"] == 0

    def test_reverse_strand(self):
        from mitoqc.modules import m5_pcg as m5
        fwd = _clean_cds(900)
        rc  = str(Seq(fwd).reverse_complement())
        seq = "A"*1000 + rc + "A"*15000
        feat = _cds("cox1", 1000, 1000+len(rc), strand=-1)
        rec = {
            "accession": "TEST", "organism": "Homo sapiens",
            "lineage": "Eukaryota; Metazoa; Chordata; Mammalia",
            "taxon_group": "Mammalia", "seq_len": len(seq), "seq": seq,
            "date": "01-JAN-2024", "platform": "Illumina",
            "features": [feat], "n_features": 1,
        }
        r = m5.run(rec)
        assert r["m5_pass"] is True

    def test_aves_nad3_exception(self):
        """Avian nad3 internal stop should be exempt (RNA-editing frameshift)."""
        from mitoqc.modules import m5_pcg as m5
        cds = "ATG" + "AAA"*4 + "TAA" + "AAA"*10 + "TAA"
        seq = "A"*1000 + cds + "A"*(20000 - 1000 - len(cds))
        feat = _cds("nad3", 1000, 1000+len(cds))
        rec = {
            "accession": "AVES001", "organism": "Gallus gallus",
            "lineage": "Eukaryota; Metazoa; Chordata; Aves",
            "taxon_group": "Aves", "seq_len": len(seq), "seq": seq,
            "date": "01-JAN-2024", "platform": "Illumina",
            "features": [feat], "n_features": 1,
        }
        r = m5.run(rec)
        assert r["m5_pass"] is True


# ── Scorer tests ──────────────────────────────────────────────────────────────

class TestScorer:
    def test_gold_all_pass(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":True,"m2_pass":True,"m3_pass":True,
                   "m4_pass":True,"m5_pass":True})
        assert r["quality_grade"] == "Gold"
        assert r["quality_score"] == 1.0

    def test_silver_one_non_critical_fail(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":True,"m2_pass":False,"m3_pass":True,
                   "m4_pass":True,"m5_pass":True})
        assert r["quality_grade"] == "Silver"

    def test_bronze_one_critical_fail(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":False,"m2_pass":True,"m3_pass":True,
                   "m4_pass":True,"m5_pass":True})
        assert r["quality_grade"] == "Bronze"

    def test_bronze_two_non_critical_fails(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":True,"m2_pass":False,"m3_pass":False,
                   "m4_pass":True,"m5_pass":True})
        assert r["quality_grade"] == "Bronze"

    def test_fail_three_modules(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":False,"m2_pass":False,"m3_pass":False,
                   "m4_pass":True,"m5_pass":True})
        assert r["quality_grade"] == "Fail"

    def test_fail_both_critical(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":False,"m2_pass":True,"m3_pass":True,
                   "m4_pass":True,"m5_pass":False})
        assert r["quality_grade"] == "Bronze"

    def test_quality_score(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":True,"m2_pass":True,"m3_pass":True,
                   "m4_pass":False,"m5_pass":False})
        assert r["quality_score"] == 0.6

    def test_missing_module_not_penalised(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":True,"m2_pass":True,"m3_pass":True,
                   "m4_pass":True})  # m5_pass missing
        assert r["quality_grade"] == "Gold"

    def test_custom_critical_modules(self):
        from mitoqc.modules.scorer import grade
        r = grade({"m1_pass":True,"m2_pass":False,"m3_pass":True,
                   "m4_pass":True,"m5_pass":True},
                  config={"critical_modules": {"m2_pass"}})
        assert r["quality_grade"] == "Bronze"


# ── Taxonomy tests ────────────────────────────────────────────────────────────

class TestTaxonomy:
    def test_mammalia(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Chordata; Mammalia; Primates") == "Mammalia"

    def test_aves(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Chordata; Aves; Passeriformes") == "Aves"

    def test_reptilia_lepidosauria(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Chordata; Lepidosauria; Squamata") == "Reptilia"

    def test_insecta(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Arthropoda; Insecta; Diptera") == "Insecta"

    def test_arachnida(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Arthropoda; Arachnida; Araneae") == "Arachnida"

    def test_crustacea_malacostraca(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Arthropoda; Malacostraca; Decapoda") == "Crustacea"

    def test_nematoda(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Nematoda; Chromadorea") == "Nematoda"

    def test_echinodermata(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Echinodermata; Asteroidea") == "Echinodermata"

    def test_other_metazoa_fallback(self):
        from mitoqc.utils.taxonomy import classify_taxon
        assert classify_taxon("Eukaryota; Metazoa; Platyhelminthes; Trematoda") == "Other_Metazoa"

    def test_insecta_not_misclassified_as_arachnida(self):
        from mitoqc.utils.taxonomy import classify_taxon
        lin = "Eukaryota; Metazoa; Arthropoda; Insecta; Coleoptera"
        assert classify_taxon(lin) == "Insecta"
