FROM python:3.7.10

WORKDIR C:/Code/WP1.3_multiple_modalities/OOP/app

COPY requirements.txt ./

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD [ "python", "./app.py C:\Code\WP1.3_multiple_modalities\OOP\configs NCITA_TEST_RR2" ]