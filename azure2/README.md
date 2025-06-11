# The azure2 directory

When you run `scripts/run_gui.sh` to run the GUI in a docker container, the script mounts this directory.

Users can put .azr files directly into `azure2` and they will be available to load in the GUI.  By default, AZURE2 will expect data to be in the `azure2/data` directory and will put outputs into the `azure2/outputs` directory.

If you aren't running the docker scripts you can ignore this directory entirely.
