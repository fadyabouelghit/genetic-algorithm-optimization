function notify_results(subject, body)
%NOTIFY_RESULTS  Dispatch a result block to configured channels.
%
% All destinations and credentials are read from outside this repo:
%   * `.notify_local.json`  (gitignored, sits next to this file) holds the
%                            ntfy topic and email addresses/server.
%   * MATLAB user prefs     hold the SMTP username/password for email
%                           (set once via `setpref('Internet', ...)`).
%
% The repo only ships `.notify_local.example.json` -- copy it to
% `.notify_local.json` and edit your settings there. If no config file
% exists, this function is a silent no-op so the rest of your code never
% breaks just because notifications aren't set up yet.
%
% Failures (network, auth, bad config) are logged via warning() but never
% thrown -- a failed ping shouldn't tank a 16-minute GA run.
%
% USAGE
%   notify_results('GA run complete', summaryString);
%
% .notify_local.json schema (see .notify_local.example.json):
%   {
%     "ntfy":  { "enabled": true,
%                "topic":   "<random-string-only-you-know>" },
%     "email": { "enabled":     false,
%                "to":          "you@example.com",
%                "from":        "ga-notifier@example.com",
%                "smtpServer":  "smtp.gmail.com",
%                "tls":         true,
%                "port":        465 }
%   }
%
% Email password setup (run ONCE outside this repo, e.g. in a startup
% script in your home directory or interactively in the MATLAB command
% window):
%   setpref('Internet', 'SMTP_Username', 'you@gmail.com');
%   setpref('Internet', 'SMTP_Password', 'your-app-password');
% Use a Google "app password" for Gmail, not your account password.

    cfg = local_load_config();
    if isempty(cfg)
        return;   % no config -> intentional silent no-op
    end

    if local_enabled(cfg, 'ntfy')
        try
            local_send_ntfy(cfg.ntfy, subject, body);
        catch ME
            warning('notify_results:ntfy', 'ntfy delivery failed: %s', ME.message);
        end
    end

    if local_enabled(cfg, 'email')
        try
            local_send_email(cfg.email, subject, body);
        catch ME
            warning('notify_results:email', 'email delivery failed: %s', ME.message);
        end
    end
end


function cfg = local_load_config()
    here = fileparts(mfilename('fullpath'));
    pth = fullfile(here, '.notify_local.json');
    if ~exist(pth, 'file')
        cfg = [];
        return;
    end
    try
        cfg = jsondecode(fileread(pth));
    catch ME
        warning('notify_results:badConfig', ...
            '.notify_local.json could not be parsed: %s', ME.message);
        cfg = [];
    end
end


function tf = local_enabled(cfg, name)
    tf = isfield(cfg, name) && isstruct(cfg.(name)) ...
        && isfield(cfg.(name), 'enabled') && logical(cfg.(name).enabled);
end


function local_send_ntfy(opts, subject, body)
    if ~isfield(opts, 'topic') || isempty(opts.topic)
        error('notify_results:ntfy:noTopic', 'ntfy.topic is empty.');
    end
    url = sprintf('https://ntfy.sh/%s', char(opts.topic));
    options = weboptions( ...
        'MediaType',         'text/plain', ...
        'CharacterEncoding', 'UTF-8', ...
        'HeaderFields',      {'Title', char(subject)});
    webwrite(url, char(body), options);
end


function local_send_email(opts, subject, body)
    required = {'to', 'from', 'smtpServer'};
    for k = 1:numel(required)
        if ~isfield(opts, required{k}) || isempty(opts.(required{k}))
            error('notify_results:email:cfg', 'email.%s missing.', required{k});
        end
    end

    setpref('Internet', 'E_mail',      char(opts.from));
    setpref('Internet', 'SMTP_Server', char(opts.smtpServer));

    if isfield(opts, 'tls') && logical(opts.tls)
        port = 465;
        if isfield(opts, 'port') && ~isempty(opts.port)
            port = opts.port;
        end
        props = java.lang.System.getProperties();
        props.setProperty('mail.smtp.auth',                'true');
        props.setProperty('mail.smtp.starttls.enable',     'true');
        props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
        props.setProperty('mail.smtp.socketFactory.port',  num2str(port));
    end

    sendmail(char(opts.to), char(subject), char(body));
end
