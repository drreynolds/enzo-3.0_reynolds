//
// Field Collection Tests
//
// Authors: Matthew Turk
//          Greg Bryan

#include "gtest/gtest.h"
#include "FieldObjects.h"

namespace {
  class GridSimpleTest : public ::testing::Test {
    protected:
      virtual void SetUp() {
        this->g = new Grid();
      }

      virtual void TearDown() {
        delete this->g;
      }

      Grid *g;

  };

}

TEST_F(GridSimpleTest, TestFieldManipulation) {
  this->g->AddField("Density", "g / cc", 2, InterpolateDirectly, CellCentered);
  this->g->AddField("H2I_Density", "g / cc", 2, InterpolateDirectly, CellCentered);

  FieldDescriptor *rho = g->GetField("Density");
  FieldDescriptor *H2 = g->GetField("H2I_Density");

  rho->CopyFrom(1e-23);
  H2->CopyFrom(1e-24);

  FieldDescriptor *H2f = H2->Duplicate();
  (*H2f) /= rho;

  ASSERT_DOUBLE_EQ(H2f->min(), 0.1);
  ASSERT_DOUBLE_EQ(H2f->max(), 0.1);
  
  // NOTE: Because we convert these from pointers to local values, they are
  // de-allocated as soon as this routine returns.

  delete H2f;
}

TEST_F(GridSimpleTest, TestFieldAccess) {
  int Dims[MAX_DIMENSIONS] = {16, 16, 16};
  long long LeftEdge[MAX_DIMENSIONS] = {0, 0, 0};
  FieldDescriptor *fd1 = new FieldDescriptor(
      CellCentered, 3, Dims, LeftEdge,
      InterpolateDirectly, "Density", "g/cc", NULL);
  fd1->CopyFrom(2.0);
  this->g->AttachField("Density", fd1);
  FieldDescriptor *fd2 = this->g->GetField("Density");
  // This tests that fd2 points to the same stuff as fd1
  ASSERT_EQ(fd2->max(), 2.0);
  ASSERT_EQ(fd2->min(), 2.0);
  // Let's see if our changes propagate...
  (*fd1) += 1.0;
  ASSERT_EQ(fd2->max(), 3.0);
  ASSERT_EQ(fd2->min(), 3.0);
}

TEST_F(GridSimpleTest, TestFieldCreation) {
  FieldDescriptor *fd = this->g->AddField("Density", "g / cc", 2, InterpolateDirectly, CellCentered);
  ASSERT_STREQ(fd->GetName(), "Density");
  ASSERT_EQ(fd->GetSize(), 20*20*20);
  long long Position[MAX_DIMENSIONS];
  fd->GetLeftEdge(Position);
  ASSERT_EQ(Position[0], -2);
  ASSERT_EQ(Position[1], -2);
  ASSERT_EQ(Position[2], -2);
  int Dimensions[MAX_DIMENSIONS];
  fd->GetCellDimensions(Dimensions);
  ASSERT_EQ(Dimensions[0], 20);
  ASSERT_EQ(Dimensions[1], 20);
  ASSERT_EQ(Dimensions[2], 20);
}

TEST(GridBaryonFieldsTest, TestFieldCreation) {
  Grid *g = new Grid();
  g->AddField("Something1", "", 0, MultiplyByDensity, CellCentered);
  FieldDescriptor *fd = g->GetField("Something1");
  ASSERT_EQ(fd->GetValues(), g->BaryonFields[0]);
  double **v = g->BaryonFields;
  delete fd;
  ASSERT_TRUE(g->BaryonFields[0] == NULL);
  g->Fields["Something1"] = NULL;
  // Note that at this time, our Grid is in a bad state!
  delete g;
}

TEST(GridBaryonFieldTest, TestFieldManagement) {
  // This is just a valgrind test
  Grid *g = new Grid();
  g->AddField("Something1", "", 0, MultiplyByDensity, CellCentered);
  double **v = g->BaryonFields;
  delete g;
}

TEST(GridBaryonFieldTest, TestAllocateOutFromUnderYou) {
  Grid *g = new Grid();
  double *v;
  g->AddField("Something1", "", 0, MultiplyByDensity, CellCentered);
  g->AddField("Something2", "", 0, MultiplyByDensity, CellCentered);
  delete g->BaryonFields[0];
  g->BaryonFields[0] = NULL;
  ASSERT_TRUE(g->GetField("Something1")->GetValues() == NULL);
  ASSERT_TRUE(g->GetField("Something2")->GetValues() == g->BaryonFields[1]);
  v = g->BaryonFields[0] = new double[16*16*16];
  ASSERT_TRUE(g->GetField("Something1")->GetValues() == v);
  ASSERT_TRUE(g->GetField("Something2")->GetValues() == g->BaryonFields[1]);
  delete g;
}

TEST(GridBaryonFieldTest, TestFindField) {
  Grid *g = new Grid();
  
  delete g;
}
