from pytest import mark


@mark.xfail
def test_run():
    assert 0


@mark.xfail
def test_get_homo_lumo_gap():
    assert 0
