import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import pytest


# To test:
# python -m pytest tests/test_all.py -v
#
# To get code coverage:
# python -m pytest --cov=coordinate_matching tests/test_all.py
#
# To get code coverage html report to see which lines we missed:
# python -m pytest --cov=coordinate_matching tests/test_all.py --cov-report=html




class TestMatchCatalogs:

    def test_simple_example(self):
        """Test a simple example for matching two catalogs."""
        from coordinate_matching import match_catalogs

        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2) # limit can vary, depending on how conservative you want to be.

        expected_idx = [0, 1]


        assert list(idx) == expected_idx

    def test_simple_example_reversed(self):
        """Test a simple example for matching two catalogs."""
        from coordinate_matching import match_catalogs

        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [2+epsilon,1-epsilon], [5+epsilon,4-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2) # limit can vary, depending on how conservative you want to be.

        expected_idx = [1, 0]


        assert list(idx) == expected_idx


    def test_simple_example_out_of_bounds(self):
            """Test a simple example for matching two catalogs where one element is too far."""
            from coordinate_matching import match_catalogs
    
            epsilon = 1e-4
            ra1, dec1 = [1,2,3], [4,5,6]
            ra2, dec2 = [1+epsilon,2-epsilon, 10], [4+epsilon,5-epsilon,20]
    
            df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
            df2 = pd.DataFrame({'ra':ra2,'dec':dec2})
    
            idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2) # limit can vary, depending on how conservative you want to be.
    
            expected_idx = [0, 1, np.nan]
    
            assert np.allclose(idx, expected_idx, equal_nan=True)



    def test_simple_example_remove_duplicates_false(self):
        """Test a simple example for matching two catalogs focusing on duplicates."""
        from coordinate_matching import match_catalogs

        epsilon = 1e-4
        ra1, dec1 = [1+epsilon/2,1+epsilon], [4+epsilon/2,4+epsilon]
        ra2, dec2 = [1,2], [4,5]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2, remove_duplicates = False) # limit can vary, depending on how conservative you want to be.

        expected_idx = [0, 0]

        assert np.allclose(idx, expected_idx, equal_nan=True)


    def test_simple_example_remove_duplicates_true(self):
        """Test a simple example for matching two catalogs focusing on duplicates."""
        from coordinate_matching import match_catalogs

        epsilon = 1e-4
        ra1, dec1 = [1+epsilon/2,1+epsilon], [4+epsilon/2,4+epsilon]
        ra2, dec2 = [1,2], [4,5]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2, remove_duplicates = True) # limit can vary, depending on how conservative you want to be.

        expected_idx = [0, np.nan]

        assert np.allclose(idx, expected_idx, equal_nan=True)



    def test_simple_example_recursive(self):
        """Test a simple example for matching two catalogs focusing on the recursive argument."""
        from coordinate_matching import match_catalogs

        limit = 2
        epsilon = limit/60/60

        ra1, dec1 = [1+epsilon/4,1+epsilon/3], [4+epsilon/4,4+epsilon/3]
        ra2, dec2 = [1,1+epsilon], [4,4+epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = limit, remove_duplicates = True, recursive = True) # limit can vary, depending on how conservative you want to be.

        expected_idx = [0, 1]
        assert np.allclose(idx, expected_idx, equal_nan=True)
        



    def test_simple_example_verbose(self):
        """Test a simple example for matching two catalogs for the verbose argument."""
        from coordinate_matching import match_catalogs

        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2, verbose = True) # limit can vary, depending on how conservative you want to be.

        expected_idx = [0, 1]


        assert list(idx) == expected_idx


class TestMergeCatalogs:

    def test_simple_example(self):
        """Test a simple example for merging two catalogs."""
        from coordinate_matching import match_catalogs, merge_catalogs


        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        idx, _, _ = match_catalogs([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], limit = 2) # limit can vary, depending on how conservative you want to be.


        df_merged = merge_catalogs(df1,df2,idx)

        expected_df_merged = pd.DataFrame({'ra_A':ra1,'dec_A':dec1,'ra_B':ra2,'dec_B':dec2})

        assert df_merged.equals(expected_df_merged)













class TestSexagesimalToDegree:
    

    def test_valid_sexagesimal_conversion(self):
        """Test converting standard RA (hms) and Dec (dms) strings."""
        from coordinate_matching import sexagesimal_to_degree

        ra_str = "12h30m00s"
        dec_str = "+45d00m00s"
        
        # Expected: 12.5 hours = 187.5 degrees, 45 degrees = 45.0 degrees
        expected_ra = 187.5
        expected_dec = 45.0
        
        ra_deg, dec_deg = sexagesimal_to_degree(ra_str, dec_str)
        
        assert ra_deg == pytest.approx(expected_ra)
        assert dec_deg == pytest.approx(expected_dec)

    def test_known_astrophysical_object(self):
        """Test with known coordinates (e.g., Crab Nebula: 05h34m31.94s, +22d00m52.2s)."""
        from coordinate_matching import sexagesimal_to_degree

        ra_str = "05h34m31.94s"
        dec_str = "+22d00m52.2s"
        
        ra_deg, dec_deg = sexagesimal_to_degree(ra_str, dec_str)
        
        # Verify result matches direct Astropy SkyCoord calculation
        expected_coord = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))
        
        assert ra_deg == pytest.approx(expected_coord.ra.degree)
        assert dec_deg == pytest.approx(expected_coord.dec.degree)

    def test_negative_declination(self):
        """Test that negative declinations are parsed correctly."""
        from coordinate_matching import sexagesimal_to_degree

        ra_str = "00h00m00s"
        dec_str = "-30d15m00s"
        
        # Expected Dec: -30 - (15/60) = -30.25 degrees
        ra_deg, dec_deg = sexagesimal_to_degree(ra_str, dec_str)
        
        assert ra_deg == pytest.approx(0.0)
        assert dec_deg == pytest.approx(-30.25)






class TestPlottingFunctions:
    """
    Not comparing/testing the actual plot, but making sure the code actually runs.
    """
    @pytest.fixture(autouse=True)
    def cleanup_plots(self, monkeypatch):
        """Automatically close all figures after each test to prevent memory leaks."""

        matplotlib.use('Agg')
        monkeypatch.setattr(plt, "show", lambda: None)

        yield

        plt.close('all')
        

        

    def test_plot_separation_limit(self):
        """Test the plot_separation_limit function."""
        from coordinate_matching import plot_separation_limit

        
        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        plot_separation_limit([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], logscale = False)



    def test_plot_separation_limit_logscale(self):
        """Test the plot_separation_limit function."""
        from coordinate_matching import plot_separation_limit

        
        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        plot_separation_limit([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], logscale = True)





    def test_plot_on_sky(self):
        """Test the plot_on_sky function."""
        from coordinate_matching import plot_on_sky 

        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        plot_on_sky([df1['ra'],df1['dec']],[df2['ra'],df2['dec']])


    
    def test_plot_on_sky_fracplot(self):
        """Test the plot_on_sky function."""
        from coordinate_matching import plot_on_sky 

        epsilon = 1e-4
        ra1, dec1 = [1,2,3], [4,5,6]
        ra2, dec2 = [1+epsilon,2-epsilon,3], [4+epsilon,5-epsilon,6]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        plot_on_sky([df1['ra'],df1['dec']],[df2['ra'],df2['dec']], frac_plot = 0.8)


    def test_plot_coordinate_difference(self):
        """Test the plot_coordinate_difference function."""
        from coordinate_matching import plot_coordinate_difference  

        epsilon = 1e-4
        ra1, dec1 = [1,2], [4,5]
        ra2, dec2 = [1+epsilon,2-epsilon], [4+epsilon,5-epsilon]

        df1 = pd.DataFrame({'ra':ra1,'dec':dec1})
        df2 = pd.DataFrame({'ra':ra2,'dec':dec2})

        plot_coordinate_difference([df1['ra'],df1['dec']],[df2['ra'],df2['dec']])


        