#include "gtest/gtest.h"
#include <hemocell.h>

// illustrative test: assert hemocell exits with non-existing error code
TEST(HemoCell, init)
{
        char *args[] = {(char*)"test", (char*)"path", NULL};
        char *inp = (char*)"config.xml";
        EXPECT_EXIT(hemo::HemoCell hemocell(inp, 0, args), ::testing::ExitedWithCode(1), "");

}
