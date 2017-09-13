import sys
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

def send_email_when_finished(set_, thoipapyoutput_file_loc,output_png_loc):
    """ Sends an email to specified address when job is finished

    Parameters
    ----------
    s : dict
        Settings dictionary extracted from excel settings file.

    Returns
    -------
    nothing but sends an email

    """

    #fromaddr = "***REMOVED***"
    fromaddr = "***REMOVED***"
    #toaddr = "***REMOVED***"
    toaddr = set_["email_to"]

    msg = MIMEMultipart()

    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = "thoipapy prediction is finished"
    body = '{a}\n\n processed list: {b}'.format(a='Dear users:', b='Your submitted job was finished, the prediciton result and fitted sine curve are in the attachment files')
    msg.attach(MIMEText(body, 'plain'))

    email_fig_list = []
    # for settings_parameter in s.keys():
    #     if settings_parameter[-5:] == "email":
    #         if s[settings_parameter] == True:
    #             email_fig_list.append(settings_parameter)

    # if email_fig_list != []:
    #     for email_fig in email_fig_list:
    #         Fig_name = email_fig[:-6]
    #         filepath = os.path.join(pathdict["single_list_fig_path"], Fig_name + ".png")
    #         print("filepath", filepath)
    #         if os.path.isfile(filepath):
    #             attachment = open(filepath, "rb")
    #             part = MIMEBase('application', 'octet-stream')
    #             part.set_payload((attachment).read())
    #             encoders.encode_base64(part)
    #             part.add_header('Content-Disposition', "attachment; filename= %s" % Fig_name + ".png")
    #             msg.attach(part)

    #filepath="/home/students/zeng/workspace/test2/out/58a2cc9b7ae16/output.csv"
    filepath=output_file_loc
    pngpath=output_png_loc
    attachment_csv = open(filepath, "rb")
    
    part = MIMEBase('application', 'octet-stream')
    part.set_payload((attachment_csv).read())
    
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', "attachment_csv; filename= output.csv ")

    msg.attach(part)
    attachment_png=open(pngpath,"rb")
    part = MIMEBase('application', 'octet-stream')
    part.set_payload((attachment_png).read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', "attachment_png; filename= output.png ")
    msg.attach(part)

    server = smtplib.SMTP('smtp.gmail.com', 587)
    #server = smtplib.SMTP('smtp.mail.com',587)
    server.starttls()
    server.login("***REMOVED***", "***REMOVED***")
    #server.login("***REMOVED***", "***REMOVED***")

    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()
    sys.stdout.write('Email sent to {}'.format(toaddr))
    attachment_csv.close()
    attachment_png.close()