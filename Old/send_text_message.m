function send_text_message(subject,message)
%    send_text_message('subject','message') or send_text_message('message')

% =========================================================================
mail = 'emptyeggcartons@gmail.com';    %GMail email address
password = 'peaches_4542';          %GMail password
% =========================================================================

if nargin == 1
    message = subject;
    subject = '';
end

emailto = '7013308270@vtext.com';

%% Set up Gmail SMTP service.
% Set up preferences
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%% Send the email
sendmail(emailto,subject,message)
