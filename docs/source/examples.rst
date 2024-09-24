Examples
=============

Installation/Usage:
*******************

::

    pip install rda-toolbox

or

::

    pip install "git+https://github.com/Robotic-Discovery-of-Antiinfectives/rda-toolbox.git"


Read raw data:
*******************

.. code-block:: Python
    :linenos:

    #!/usr/bin/env python3

    import rda_toolbox as rda
    import glob

    rda.readerfiles_rawdf(glob.glob("path/to/raw/readerfiles/*"))


If you have multiple timepoints:

.. code-block:: Python
    :linenos:

    import rda_toolbox as rda
    import os

    timepoints_rawdata = dict()
    # change "raw/data/path" to your raw data location
    for timepoint in [f.path for f in os.scandir("raw/data/path") if f.is_dir()]:
        timepoints_rawdata[os.path.basename(timepoint)] = rda.readerfiles_rawdf(
            [os.path.join(timepoint, file) for file in os.listdir(timepoint)]
        )

Visualize:
*******************
