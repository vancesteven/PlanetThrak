function ampersand_delimmed_numstr = num2tex(num,format,varargin)
    t_dec = num2str(mod(abs(num),1),format);
    t_num = num2str(fix(num));
if strcmp(varargin,'bold')
    ampersand_delimmed_numstr = ['\textbf{' t_num  '} & \textbf{' t_dec(3:end) '}'];
else
    ampersand_delimmed_numstr = [t_num '&' t_dec(3:end)];
end
