/*************************************************************************
 ** Sample3
 ** Load sample.cfg and try to add a setting "foo"..
 **   on success save to testfoo.cfg
 *************************************************************************/

#include <iostream>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;

/***************************************************************************/

int main()
{
  Config cfg;
  try
  {
    /* Load the configuration.. */
    cout << "loading [sample.cfg]..";
    cfg.readFile("sample.cfg");
    cout << "ok" << endl;

    /* Add setting "foo" */
    cout << "add setting \"foo\"..";
    Setting &root = cfg.getRoot();
    Setting &foo  = root.add("foo",TypeInt);
    foo = 1234;
    cout << "ok" << endl;

    /* Save to "samplefoo.cfg" */
    cout << "saving [samplefoo.cfg]..";
    cfg.writeFile("samplefoo.cfg");
    cout << "ok" << endl;

    cout << "Done!" << endl;
  }
  catch (...)
  {
    cout << "failed" << endl;
  }

  return 0;
}

/***************************************************************************/
