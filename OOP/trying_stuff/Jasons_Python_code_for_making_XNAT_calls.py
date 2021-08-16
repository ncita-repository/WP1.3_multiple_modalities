""" Jason's Python code for making XNAT REST calls. 

Provided via Slack on 02/03/2021.


So the bottom part is a context manager (with statement) and means the 'xnat_session' variable keeps a live connection as a single session you can use to do anything on the xapi with
you can give the username/password combo once at the start with .auth
then using the .get, .post, .delete, etc... requests functions on the session you can give a URL http call to wherever you need

in the code above to get a list of scans IDs for a session I make the request, get the json() from it and then search for the information I care about which is in ResultSet --> Result --> ID  (there is probably a cleaner way to do this part but I haven't spent anytime thinking about it, it gets the job done anyway)

Then I download the files for a scan resource as a zip using the optional ?format=zip and save the .content of the request anywhere on the local disk using the .write (you can use the zipfile package to unzip this with python as well and save it somewhere else to work with the data using python)

https://wiki.xnat.org/display/XAPI/XNAT+REST+API+Directory

This wiki catalogues the calls you can make in decent detail, you will probably still need a bit of trial and error to make sure a call is doing what you want it to
"""


import requests


def func1():
    # do stuff
	
	
def func2():
    # more stuff
	
	
def main():
    # more code
    func1()
	
    func2()
	
    scan_id_list = []
    scans_request = xnat_session.get(f"{self.domain}/data/experiments/{self.experiment_id}/scans")
    scan_list_json = scans_request.json()
    [scan_id_list.append(scan['ID']) for scan in scan_list_json['ResultSet']['Result']]
    
	scan_full_download = xnat_session.get(f"{self.domain}/data/experiments/{self.experiment_id}/scans/{scan}/"
                                    f"resources/DICOM/files?format=zip")
									
    with open("path/to/save/data/at/", 'wb') as file:
        file.write(scan_full_download.content)
		
		
with requests.Session() as xnat_session:
    xnat_session.auth = (f'{username}', f'{password}')
    main()