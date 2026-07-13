// Fast, standalone unit tests for the pure lambda-weighting math used by
// all backmap styles (backmap_lambda_weights.h). No LAMMPS headers, no
// LAMMPS build required — this is the pre-simulation correctness gate.

#include <gtest/gtest.h>

#include "backmap_lambda_weights.h"

using LAMMPS_NS::BackmapLambda::clamp_lambda;
using LAMMPS_NS::BackmapLambda::compute_weight3;
using LAMMPS_NS::BackmapLambda::is_almost_zero;
using LAMMPS_NS::BackmapLambda::same_bead;

namespace {

// ---- CG weight: w = 1 - lambda_global, independent of same_bead ----

TEST(ComputeWeight3, CgAtLambdaZeroIsFullStrength) {
  EXPECT_DOUBLE_EQ(compute_weight3(/*same_bead=*/false, /*is_cg=*/true, 0.0),
                   1.0);
  EXPECT_DOUBLE_EQ(compute_weight3(/*same_bead=*/true, /*is_cg=*/true, 0.0),
                   1.0);
}

TEST(ComputeWeight3, CgAtLambdaOneIsZero) {
  EXPECT_DOUBLE_EQ(compute_weight3(false, true, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(compute_weight3(true, true, 1.0), 0.0);
}

TEST(ComputeWeight3, CgIsLinearInLambda) {
  EXPECT_DOUBLE_EQ(compute_weight3(false, true, 0.25), 0.75);
  EXPECT_DOUBLE_EQ(compute_weight3(false, true, 0.5), 0.5);
  EXPECT_DOUBLE_EQ(compute_weight3(false, true, 0.75), 0.25);
}

TEST(ComputeWeight3, CgPhase1OverridesToFullStrength) {
  EXPECT_DOUBLE_EQ(compute_weight3(false, true, 0.5, /*phase=*/1), 1.0);
  EXPECT_DOUBLE_EQ(compute_weight3(false, true, 1.0, /*phase=*/1), 1.0);
}

// ---- AT intra-bead (same_bead=true): always full strength, independent
// of lambda_global. Represents real intra-molecular chemistry that exists
// regardless of the CG/AT resolution ramp (matches the original
// ESPResSo++ production driver, which builds the full AT interaction list
// once, unconditionally, before the resolution ramp is ever activated). ----

TEST(ComputeWeight3, AtIntraBeadIsFullStrengthAtLambdaZero) {
  EXPECT_DOUBLE_EQ(compute_weight3(/*same_bead=*/true, /*is_cg=*/false, 0.0),
                   1.0);
}

TEST(ComputeWeight3, AtIntraBeadIsFullStrengthAtAnyLambda) {
  EXPECT_DOUBLE_EQ(compute_weight3(true, false, 1.0e-9), 1.0);
  EXPECT_DOUBLE_EQ(compute_weight3(true, false, 0.5), 1.0);
  EXPECT_DOUBLE_EQ(compute_weight3(true, false, 1.0), 1.0);
}

// ---- AT inter-bead (same_bead=false): w = lambda_global exactly ----

TEST(ComputeWeight3, AtInterBeadIsExactlyLambdaGlobal) {
  EXPECT_DOUBLE_EQ(compute_weight3(/*same_bead=*/false, /*is_cg=*/false, 0.0),
                   0.0);
  EXPECT_DOUBLE_EQ(compute_weight3(false, false, 0.3), 0.3);
  EXPECT_DOUBLE_EQ(compute_weight3(false, false, 0.5), 0.5);
  EXPECT_DOUBLE_EQ(compute_weight3(false, false, 1.0), 1.0);
}

// ---- same_bead(): pair/angle/dihedral overloads ----

TEST(SameBead, PairTrueWhenBothMapToSameBead) {
  int atom2cg[] = {5, 5, 7};
  EXPECT_TRUE(same_bead(atom2cg, 0, 1));
}

TEST(SameBead, PairFalseWhenBeadsDiffer) {
  int atom2cg[] = {5, 7};
  EXPECT_FALSE(same_bead(atom2cg, 0, 1));
}

TEST(SameBead, PairFalseWhenEitherAtomUnmapped) {
  int atom2cg[] = {-1, 5};
  EXPECT_FALSE(same_bead(atom2cg, 0, 1));
  int atom2cg2[] = {5, -1};
  EXPECT_FALSE(same_bead(atom2cg2, 0, 1));
}

TEST(SameBead, AngleTrueWhenAllThreeMapToSameBead) {
  int atom2cg[] = {3, 3, 3};
  EXPECT_TRUE(same_bead(atom2cg, 0, 1, 2));
}

TEST(SameBead, AngleFalseWhenAllThreeDiffer) {
  int atom2cg[] = {1, 2, 3};
  EXPECT_FALSE(same_bead(atom2cg, 0, 1, 2));
}

TEST(SameBead, AngleFalseWhenMixedNotAllSameNotAllDifferent) {
  // i1, i2 same bead, i3 different — must not be misclassified as "same".
  int atom2cg[] = {4, 4, 9};
  EXPECT_FALSE(same_bead(atom2cg, 0, 1, 2));
}

TEST(SameBead, DihedralTrueWhenAllFourMapToSameBead) {
  int atom2cg[] = {2, 2, 2, 2};
  EXPECT_TRUE(same_bead(atom2cg, 0, 1, 2, 3));
}

TEST(SameBead, DihedralFalseWhenOnlyThreeOfFourMatch) {
  int atom2cg[] = {2, 2, 2, 6};
  EXPECT_FALSE(same_bead(atom2cg, 0, 1, 2, 3));
}

TEST(SameBead, DihedralFalseWhenAllFourDiffer) {
  int atom2cg[] = {1, 2, 3, 4};
  EXPECT_FALSE(same_bead(atom2cg, 0, 1, 2, 3));
}

// ---- clamp_lambda(): out-of-range values are clamped before use ----

TEST(ClampLambda, NegativeClampsToZero) {
  EXPECT_DOUBLE_EQ(clamp_lambda(-0.5), 0.0);
  EXPECT_DOUBLE_EQ(clamp_lambda(-1.0e6), 0.0);
}

TEST(ClampLambda, AboveOneClampsToOne) {
  EXPECT_DOUBLE_EQ(clamp_lambda(1.5), 1.0);
  EXPECT_DOUBLE_EQ(clamp_lambda(1.0e6), 1.0);
}

TEST(ClampLambda, InRangePassesThrough) {
  EXPECT_DOUBLE_EQ(clamp_lambda(0.0), 0.0);
  EXPECT_DOUBLE_EQ(clamp_lambda(0.42), 0.42);
  EXPECT_DOUBLE_EQ(clamp_lambda(1.0), 1.0);
}

// ---- is_almost_zero(): negligible-weight skip check ----

TEST(IsAlmostZero, TrueForExactZeroAndBelowThreshold) {
  EXPECT_TRUE(is_almost_zero(0.0));
  EXPECT_TRUE(is_almost_zero(1.0e-12));
  EXPECT_TRUE(is_almost_zero(-1.0e-12));
}

TEST(IsAlmostZero, FalseAboveThreshold) {
  EXPECT_FALSE(is_almost_zero(1.0e-6));
  EXPECT_FALSE(is_almost_zero(0.5));
}

}  // namespace
