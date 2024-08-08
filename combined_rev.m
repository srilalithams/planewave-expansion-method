desired_bandgap_min = input('Enter the desired minimum bandgap: ');
desired_bandgap_max = input('Enter the desired maximum bandgap: ');
mode = input('Enter mode (''e'' for TE mode, ''m'' for TM mode): ', 's');

if mode == 'm'
    TM_mode(desired_bandgap_min,desired_bandgap_max);

elseif mode == 'e'
    TE_mode(desired_bandgap_min,desired_bandgap_max);


else
    disp('Invalid mode selected.');
end
