
from sqtl_pipeline.analysis import get_variant_ids

def test_get_variant_ids():
    splice_graph_file = "path/to/splice_graph_file"
    junction_ids = ["chr1:12345:23456:clu_1_+_1_2_d", "chr1:34567:45678:clu_2_-_3_4_a"]
    variant_ids = get_variant_ids(splice_graph_file, junction_ids)
    assert len(variant_ids) == 2
    assert variant_ids[0] == "variant_1"
    assert variant_ids[1] == "variant_2"