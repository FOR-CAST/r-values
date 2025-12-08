# Alberta mountain pine beetle r values

## Background

Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Quisque maximus enim sit amet neque auctor pretium at sed mauris.
Duis a maximus lacus.
Sed pharetra urna euismod interdum suscipit.
Morbi tempus blandit justo, dictum dapibus dui tristique et.
Etiam id orci dictum velit ultricies viverra.
Ut elit velit, fermentum a hendrerit vitae, dictum id dui.
Nulla mollis semper eros, eget luctus leo.
Fusce euismod est at eros dignissim accumsan.
Sed a laoreet nibh.
Nam pulvinar orci sed enim fringilla pellentesque.
Nulla purus ex, rutrum sit amet elit nec, blandit posuere mi.
Curabitur dignissim nisi nunc, in porta lacus porta et.
Suspendisse potenti.
In hac habitasse platea dictumst.
Nulla enim erat, egestas vel lectus ac, facilisis ornare lorem.
Vivamus egestas turpis sit amet neque imperdiet iaculis.

Duis sollicitudin massa in aliquet mollis.
Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Nunc facilisis diam augue, nec pulvinar lorem cursus ac.
Pellentesque ac massa dolor.
Praesent egestas porttitor dictum.
Vestibulum enim dolor, fermentum ac consequat id, ornare vitae risus.
Ut semper, massa at sodales facilisis, sem diam efficitur lacus, eget aliquam lorem felis sed metus.

## Data availability

Donec quis sem in ex ullamcorper lobortis at eget sapien.
Maecenas vel accumsan orci. 
Aenean blandit ut sem nec fringilla.
Mauris luctus tortor in mauris fringilla faucibus.
Curabitur et venenatis massa, nec lobortis metus.
Nam non nunc eget mi ullamcorper feugiat in ut tortor.
Fusce porttitor leo quis tristique laoreet.
Donec scelerisque lacus a accumsan suscipit.

## Initial project setup

0. Prerequisites:

    a.  install a suitable Java Runtime Environment following instructions at
        <https://github.com/CWFC-CCFB/J4R/wiki#requirements>

1. Clone the repository: 

    ```bash
    git clone https://github.com/FOR-CAST/r-values
    ```

2. Install required packages to project library:

    ```r
    renv::restore()
    ```

3. Install PhantomJS:

    ```r
    webshot::install_phantomjs()
    ```

## Running the code

Open the primary script (`00-main.R`) and run the lines sequentially.

- select the spatial extent of the analyses by setting `run_for = "AB"` for provincial scale,
  or `run_for = "NP"` for running the national parks analyses;
- setting `extract_mdb = TRUE` requires running on Windows to be able to open the Access database files;
- setting `plot_all = TRUE` will build all plots, including intermediate ones;
- setting `rerun_all = TRUE` will rebuild all local intermediate datasets (use carefully!);
- running `01-download-data.R` will prompt for Google credentials to access several data sets;
- mdb extraction is skipped gets skipped if `extract_mdb = FALSE`;
- most of the mdb import gets skipped unless the outputs it creates are missing or `rerun_all = TRUE`;
- most of the mdb import gets skipped unless the outputs it creates are missing or `rerun_all = TRUE`;
- sourcing `03-Jasper-analyses.R` can be done after the '01-download-data.R' script has been run
  (it does not need the 01a,b nor 02 scripts to have been run);

![](workflow.png)

# References

Cooke, B. J., Brett, R., Olesinski, J., Chubaty, A. M., and A. L. Carroll (*submitted*) Weather, climate, and the rise and fall of an unprecedented outbreak of mountain pine beetle in Alberta’s Mountain Parks, 1999-2023.
