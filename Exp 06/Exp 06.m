    clear
    close all

    %% Generate PN Sequence

	samples_per_frame = 63;

	pn_seq_generator = comm.PNSequence('Polynomial',[1 0 0 0 0 1 1], ...
	    'InitialConditions',[1 1 1 1 1 1],'SamplesPerFrame',samples_per_frame);
    pn_seq = pn_seq_generator();

    Symbol_Length = 1000;
    BPSK_Mod_pn_seq = zeros([Symbol_Length samples_per_frame]);
    BPSK_Mod_pn_seq(1, :) = pn_seq * 2 - 1;
    BPSK_Mod_pn_seq = BPSK_Mod_pn_seq(:);

    t = (0 : Symbol_Length * samples_per_frame - 1)' - Symbol_Length * samples_per_frame / 2;

    figure(1)
    plot(t, BPSK_Mod_pn_seq)
    title('BPSK PN Sequence')
    %% Generate Time Domain Filters

    Filter_Length = 1 * Symbol_Length;
    beta = 1/2;

    % Case 1: Rectangle
    rect_filter = abs(t) < Filter_Length / 2;
    % Case 2: Triangle
    triangle_filter = 1 / Filter_Length * (Filter_Length - abs(t)) .* (abs(t) <= Filter_Length);
    % Case 3: Raised Cosine Pulse
    raised_cosine_pulse = zeros(size(BPSK_Mod_pn_seq));
    for i = 1:Symbol_Length * samples_per_frame
	        if abs(t(i)) == Filter_Length / beta / 2
			        raised_cosine_pulse(i) = pi / 4 * sinc(t(i) / Filter_Length);
				    else
					            raised_cosine_pulse(i) = 1 * sinc(t(i) / Filter_Length) * cos(pi * beta * t(i) / Filter_Length) / (1 - (2 * beta * t(i) / Filter_Length)^2);
						        end
    end
    % Case 4: Root Raised Cosine Pulse
    root_raised_cosine_pulse = zeros(size(BPSK_Mod_pn_seq));
    for i = 1:Symbol_Length * samples_per_frame
	        if t(i) == 0
			        root_raised_cosine_pulse(i) = 1 * (1 + beta* (4 / pi -1));
				    elseif abs(t(i)) == Filter_Length / beta / 4
					            root_raised_cosine_pulse(i) = beta  / sqrt(2) * ((1 + 2 / pi) * sin(pi / 4 / beta) + (1 - 2 / pi) * cos(pi / 4 / beta));
						        else
								        root_raised_cosine_pulse(i) = 1 * (sin(pi * t(i) * (1 - beta) / Filter_Length) + 4 * beta * t(i) / Filter_Length * cos(pi * t(i) * (1 + beta) / Filter_Length)) / (pi * t(i) / Filter_Length * (1 - (4 * beta * t(i) / Filter_Length)^2));
									    end
    end

    figure(2)
    title('Time Domain Filters')
    hold on
    plot(t, rect_filter)
    plot(t, triangle_filter)
    plot(t, raised_cosine_pulse)
    plot(t, root_raised_cosine_pulse)
    hold off
    legend('rectangle', 'triangle','raised cosine', 'root raised cosine')
    %% Applying The Filters On PN Sequence

    rect_pn = 1e10 * real(ifft(fft(BPSK_Mod_pn_seq).* fft(rect_filter))) / (Symbol_Length * samples_per_frame)^2;
    triangle_pn = 1e10 * real(ifft(fft(BPSK_Mod_pn_seq).* fft(triangle_filter))) / (Symbol_Length * samples_per_frame)^2;
    raised_cosine_pn = 1e10 * real(ifft(fft(BPSK_Mod_pn_seq).* fft(raised_cosine_pulse))) / (Symbol_Length * samples_per_frame)^2;
    root_raised_cosine_pn = 1e10 * real(ifft(fft(BPSK_Mod_pn_seq).* fft(root_raised_cosine_pulse))) / (Symbol_Length * samples_per_frame)^2;

    figure(3)
    title('Convolution With PN Sequence')
    hold on
    plot(t, rect_pn)
    plot(t, triangle_pn)
    plot(t, raised_cosine_pn)
    plot(t, root_raised_cosine_pn)
    xlim([-5e3 5e3])
    hold off
    legend('rectangle', 'triangle','raised cosine', 'root raised cosine')
    %% Eye Diagram
    
    figure(4)
    subplot(2, 2, 1)
    eye_diagram(rect_pn, 2 * Symbol_Length, Symbol_Length)
    title('Rectangle')
    subplot(2, 2, 2)
    eye_diagram(triangle_pn, 2 * Symbol_Length, Symbol_Length)
    title('Triangle')
    subplot(2, 2, 3)
    eye_diagram(raised_cosine_pn, 2 * Symbol_Length, Symbol_Length)
    title('Raised Cosine')
    subplot(2, 2, 4)
    eye_diagram(root_raised_cosine_pn, 2 * Symbol_Length, Symbol_Length)
    title('Root Raised Cosine')
    sgtitle('Eye Diagram without noise')
    %% Adding Noise
    
    noise_low = 30;
    noise_med = 20;
    noise_high = 10;
    
    Power_rect = sum(rect_pn .* rect_pn) / (Symbol_Length * samples_per_frame);
    Power_triangle = sum(triangle_pn .* triangle_pn) / (Symbol_Length * samples_per_frame);
    Power_rc = sum(raised_cosine_pn .* raised_cosine_pn) / (Symbol_Length * samples_per_frame);
    Power_rrc = sum(root_raised_cosine_pn .* root_raised_cosine_pn) / (Symbol_Length * samples_per_frame);
    
    figure(5)
    subplot(2, 2, 1)
    eye_diagram(awgn(rect_pn, noise_low, Power_rect), 2 * Symbol_Length, Symbol_Length)
    title('Rectangle')
    subplot(2, 2, 2)
    eye_diagram(awgn(triangle_pn, noise_low, Power_triangle), 2 * Symbol_Length, Symbol_Length)
    title('Triangle')
    subplot(2, 2, 3)
    eye_diagram(awgn(raised_cosine_pn, noise_low, Power_rc), 2 * Symbol_Length, Symbol_Length)
    title('Raised Cosine')
    subplot(2, 2, 4)
    eye_diagram(awgn(root_raised_cosine_pn, noise_low, Power_rrc), 2 * Symbol_Length, Symbol_Length)
    title('Root Raised Cosine')
    sgtitle("Low Noise")
    
    figure(6)
    subplot(2, 2, 1)
    eye_diagram(awgn(rect_pn, noise_med, Power_rect), 2 * Symbol_Length, Symbol_Length)
    title('Rectangle')
    subplot(2, 2, 2)
    eye_diagram(awgn(triangle_pn, noise_med, Power_triangle), 2 * Symbol_Length, Symbol_Length)
    title('Triangle')
    subplot(2, 2, 3)
    eye_diagram(awgn(raised_cosine_pn, noise_med, Power_rc), 2 * Symbol_Length, Symbol_Length)
    title('Raised Cosine')
    subplot(2, 2, 4)
    eye_diagram(awgn(root_raised_cosine_pn, noise_med, Power_rrc), 2 * Symbol_Length, Symbol_Length)
    title('Root Raised Cosine')
    sgtitle("Medium Noise")
    
    figure(7)
    subplot(2, 2, 1)
    eye_diagram(awgn(rect_pn, noise_high, Power_rect), 2 * Symbol_Length, Symbol_Length)
    title('Rectangle')
    subplot(2, 2, 2)
    eye_diagram(awgn(triangle_pn, noise_high, Power_triangle), 2 * Symbol_Length, Symbol_Length)
    title('Triangle')
    subplot(2, 2, 3)
    eye_diagram(awgn(raised_cosine_pn, noise_high, Power_rc), 2 * Symbol_Length, Symbol_Length)
    title('Raised Cosine')
    subplot(2, 2, 4)
    eye_diagram(awgn(root_raised_cosine_pn, noise_high, Power_rrc), 2 * Symbol_Length, Symbol_Length)
    title('Root Raised Cosine')
    sgtitle("High Noise")
    %% Function Eye Diagram
    
    function eye_diagram(signal, Length, Step)
        Signal_Length = length(signal);
        number_of_segments = Signal_Length / Step;
        
        hold on
        for i = 0:number_of_segments - 1
            plot(signal(i * Step + 1 : min(Signal_Length, i * Step + Length)), 'b')
        end
        hold off
    end
    
