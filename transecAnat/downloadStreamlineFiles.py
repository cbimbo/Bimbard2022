import pandas as pd
dataPath = "C:/Users/Hamish/Desktop/InjectionsA1_AllenBrainAtlas/projection_search_results.csv" # has to be downloaded beforehand
df = pd.read_csv(dataPath)

from selenium import webdriver
from selenium.webdriver.common.by import By
from webdriver_manager.chrome import ChromeDriverManager

browser = webdriver.Chrome(ChromeDriverManager().install())

browser.get('https://neuroinformatics.nl/HBP/allen-connectivity-viewer/streamline-downloader.html')
id_box = browser.find_element(By.ID,'connectivityId')

for d in df.id:
    id_box.clear()
    id_box.send_keys(str(d))
    browser.find_element(By.XPATH,'//button[text()="Retrieve the streamline data"]').click()


# might have to click on "allow" for the multiple downloads to start!
# does not download them all on the first try -- maybe add a little pause