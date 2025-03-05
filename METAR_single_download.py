# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:43:25 2024

@author: Will Maslanka

Script for scraping data from Iowa State University METAR Archive
https://mesonet.agron.iastate.edu/

Downloads all METAR data for a single station, between two given datetimes.

"""
from datetime import datetime
from urllib.request import urlopen
import os, time

def make_directory(save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
def download_data(uri):
    """Fetch the data from the IEM

    The IEM download service has some protections in place to keep the number
    of inbound requests in check.  This function implements an exponential
    backoff to keep individual downloads from erroring.

    Args:
      uri (string): URL to fetch

    Returns:
      string data
    """
    attempt = 0
    while attempt < MAX_ATTEMPTS:
        try:
            data = urlopen(uri, timeout=300).read().decode("utf-8")
            if data is not None and not data.startswith("ERROR"):
                return data
        except Exception as exp:
            print(f"download_data({uri}) failed with {exp}")
            time.sleep(5)
        attempt += 1

    print("Exhausted attempts to download, returning empty data")
    return ""

#%% User Defined Variables
# timestamps in UTC to request data for
start_dt = datetime(2024,9,27,0,0) # start datetime
end_dt = datetime(2024,9,28,0,0) # end datetime
station = 'CZOL' # IACO Airport Code
# Codes can be found here: https://www.prokerala.com/travel/airports/
# Or google

save_dir = 'C:/Users/k2262276/Documents/Canada_2024/ToUpload/'
make_directory(save_dir)

# Number of attempts to download data
MAX_ATTEMPTS = 6

# HTTPS here can be problematic for installs that don't have Lets Encrypt CA
SERVICE = "http://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
service = SERVICE + "data=all&tz=Etc/UTC&format=comma&latlon=yes&"
service += start_dt.strftime("year1=%Y&month1=%m&day1=%d&")
service += end_dt.strftime("year2=%Y&month2=%m&day2=%d&")
service += f"&station={station}"

# Download the data from the url
data = download_data(service)

# Save the data
outfn = f"{save_dir}{station}_{start_dt:%Y%m%d%H%M}_{end_dt:%Y%m%d%H%M}.txt"
with open(outfn, "w", encoding="ascii") as fh:
    fh.write(data)

