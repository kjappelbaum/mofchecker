"""Test checks on the structure graph."""
from mofchecker.checks.global_structure.graphcheck import IsThreeDimensional


def test_is_three_dimensional(get_3d_structure_and_graph, get_1d_structure_and_graph):
    """Test the IsThreeDimensional check."""
    _structure, graph = get_3d_structure_and_graph
    check = IsThreeDimensional(graph)
    assert check.is_ok
    assert check.name == "3D connected structure graph."
    assert check.description == "Check if the structure graph is 3D connected."

    _structure, graph = get_1d_structure_and_graph
    check = IsThreeDimensional(graph)
    assert not check.is_ok
    assert check.name == "3D connected structure graph."
    assert check.description == "Check if the structure graph is 3D connected."
