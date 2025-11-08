#!/usr/bin/env python3
"""
Unit tests for breakpoint parsing.
"""
import pytest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from convert_svaba_xlsx_to_bedpe import parse_breakpoint


def test_parse_breakpoint_standard():
    """Test standard breakpoint format: chrom:pos-pos(strand)"""
    # Standard format
    result = parse_breakpoint("13:33777850-33777850+")
    assert result == ("chr13", 33777849, 33777850, "+")
    
    result = parse_breakpoint("1:1000000-1000000-")
    assert result == ("chr1", 999999, 1000000, "-")


def test_parse_breakpoint_chr_prefix():
    """Test breakpoints with chr prefix"""
    result = parse_breakpoint("chr13:33777850-33777850+")
    assert result == ("chr13", 33777849, 33777850, "+")
    
    result = parse_breakpoint("chrX:123456-123456-")
    assert result == ("chrX", 123455, 123456, "-")


def test_parse_breakpoint_sex_chromosomes():
    """Test X, Y chromosomes"""
    result = parse_breakpoint("X:50000000-50000000+")
    assert result == ("chrX", 49999999, 50000000, "+")
    
    result = parse_breakpoint("Y:10000000-10000000-")
    assert result == ("chrY", 9999999, 10000000, "-")


def test_parse_breakpoint_mitochondrial():
    """Test mitochondrial chromosome"""
    result = parse_breakpoint("M:1000-1000+")
    assert result == ("chrM", 999, 1000, "+")
    
    result = parse_breakpoint("MT:2000-2000-")
    assert result == ("chrMT", 1999, 2000, "-")


def test_parse_breakpoint_whitespace():
    """Test breakpoints with whitespace"""
    result = parse_breakpoint("  13:33777850-33777850+  ")
    assert result == ("chr13", 33777849, 33777850, "+")
    
    result = parse_breakpoint("chr1:1000000-1000000-")
    assert result == ("chr1", 999999, 1000000, "-")


def test_parse_breakpoint_case_insensitive():
    """Test case-insensitive parsing"""
    result = parse_breakpoint("x:50000000-50000000+")
    assert result == ("chrX", 49999999, 50000000, "+")
    
    result = parse_breakpoint("CHR13:33777850-33777850+")
    assert result == ("chr13", 33777849, 33777850, "+")


def test_parse_breakpoint_invalid():
    """Test invalid breakpoint formats"""
    assert parse_breakpoint("") is None
    assert parse_breakpoint("invalid") is None
    assert parse_breakpoint("13:33777850") is None  # Missing strand
    assert parse_breakpoint("13:33777850-33777850") is None  # Missing strand
    assert parse_breakpoint(None) is None


def test_parse_breakpoint_zero_position():
    """Test position 1 (should give start=0, end=1)"""
    result = parse_breakpoint("1:1-1+")
    assert result == ("chr1", 0, 1, "+")


def test_parse_breakpoint_point_breakend():
    """Test that point breakends have end = start + 1"""
    result = parse_breakpoint("13:33777850-33777850+")
    chrom, start, end, strand = result
    assert end == start + 1, f"Point breakend should have end=start+1, got end={end}, start={start}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

