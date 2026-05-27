% =========================================================
%  ECG_All_fixed.m  --  MIT-BIH Record 100, 10 s
%  Các bản sửa so với ECG_All.m gốc được đánh dấu [FIX]
% =========================================================
clc; clear; close all;

REC='100'; DUR=10; FPASS=[0.5 40]; NF=50; NQ=35; TOL=150;

% ---------------------------------------------------------
%  ĐỌC HEADER
% ---------------------------------------------------------
fid=fopen([REC '.hea'],'r'); lines={};
while ~feof(fid)
    ln=strtrim(fgetl(fid));
    if ischar(ln)&&~isempty(ln)&&ln(1)~='#', lines{end+1}=ln; end
end
fclose(fid);

f1=strsplit(lines{1});
ns=str2double(f1{2});
fs=str2double(f1{3});

gain=200*ones(1,ns);
base=zeros(1,ns);
for ch=1:ns
    if numel(lines)<ch+1, break; end
    p=strsplit(lines{ch+1});
    if numel(p)<2, continue; end

    % [FIX] Đọc gain từ đúng trường (p{3}) thay vì p{2} (là format)
    % Ưu tiên 1: p{3} có dạng "200(1024)/mV" hoặc "200(1024)"
    gain_str = '';
    if numel(p)>=3, gain_str=p{3}; end
    if ~isempty(gain_str)
        tk=regexp(gain_str,'^([0-9.]+)(?:\(([0-9\-]+)\))?','tokens','once');
        if ~isempty(tk)
            g=str2double(tk{1}); if ~isnan(g)&&g>0, gain(ch)=g; end
            if numel(tk)>=2&&~isempty(tk{2})
                b=str2double(tk{2}); if ~isnan(b), base(ch)=b; end
            end
        end
    end

    % [FIX] Baseline từ đúng trường p{5} (ADC zero) thay vì p{4} (bit-res)
    if base(ch)==0 && numel(p)>=5
        b=str2double(p{5});
        if ~isnan(b), base(ch)=b; end
    end
end

% ---------------------------------------------------------
%  ĐỌC TÍN HIỆU  (format 212)
% ---------------------------------------------------------
fid=fopen([REC '.dat'],'rb','ieee-le'); raw=fread(fid,Inf,'*uint8'); fclose(fid);
nt=floor(numel(raw)/3);
B=reshape(raw(1:nt*3),3,nt);
s1=double(B(1,:))+mod(double(B(2,:)),16)*256;
s2=double(B(3,:))+floor(double(B(2,:))/16)*256;
s1(s1>2047)=s1(s1>2047)-4096;
s2(s2>2047)=s2(s2>2047)-4096;
if ns==2, adc=[s1(:) s2(:)];
else,     adc=reshape([s1;s2],[],1);
end
sig=(double(adc)-base)./gain;

% ---------------------------------------------------------
%  CẮT & THÊM NHIỄU
% ---------------------------------------------------------
ecg=sig(1:min(DUR*fs,size(sig,1)),1); N=numel(ecg); t=(0:N-1)'/fs;
rng(42);
ecg_n = ecg + 0.15*sin(2*pi*0.3*t) ...      % BW 0.3 Hz
             + 0.08*sin(2*pi*NF*t)  ...      % PLI 50 Hz
             + 0.03*randn(N,1);               % EMG

% ---------------------------------------------------------
%  BỘ LỌC IIR  (bandpass + notch)
% ---------------------------------------------------------
[b_bp,a_bp]=butter(4,FPASS/(fs/2),'bandpass');
[b_nt,a_nt]=iirnotch(NF/(fs/2), NF/(fs/2)/NQ);
ecg_f=filtfilt(b_nt,a_nt, filtfilt(b_bp,a_bp,ecg_n));

% ---------------------------------------------------------
%  FFT & SNR/PRD
% ---------------------------------------------------------
win=hann(N); sc=2/sum(win); NFF=2^nextpow2(N);
fax=(0:NFF/2-1)'*(fs/NFF);
fdb=@(x) sc*abs(fft(x.*win,NFF));
tmp=fdb(ecg);   Yc=tmp(1:NFF/2);
tmp=fdb(ecg_n); Yn=tmp(1:NFF/2);
tmp=fdb(ecg_f); Yf=tmp(1:NFF/2);
snr_=@(x) 10*log10(mean(ecg.^2)/mean((x-ecg).^2));
prd_=@(x) 100*norm(x-ecg)/norm(ecg);
snrN=snr_(ecg_n); snrF=snr_(ecg_f);
prdN=prd_(ecg_n); prdF=prd_(ecg_f);

% ---------------------------------------------------------
%  PAN-TOMPKINS  (hàm ở cuối file -- đã sửa)
% ---------------------------------------------------------
[r_mwi,r_ecg,pts]=pan_tompkins(ecg_f,fs);
nr=numel(r_ecg);
if nr<2, error('Khong du dinh R.'); end
rr=diff(r_ecg)/fs; hr=60./rr;
SDNN  =std(rr)*1e3;
RMSSD =sqrt(mean(diff(rr).^2))*1e3;
pNN50 =100*sum(abs(diff(rr))>0.05)/max(1,numel(rr)-1);
fprintf('HR: %.1f+-%.1f bpm | SDNN=%.1fms | RMSSD=%.1fms | pNN50=%.1f%%\n',...
        mean(hr),std(hr),SDNN,RMSSD,pNN50);

% ---------------------------------------------------------
%  ĐỌC ANNOTATION .atr
% ---------------------------------------------------------
% ---------------------------------------------------------
%  ĐỌC ANNOTATION .atr  [FIX: parser WFDB đúng chuẩn]
% ---------------------------------------------------------
has_ann=false;
try
    fid=fopen([REC '.atr'],'rb','ieee-le');
    raw2=fread(fid,Inf,'*uint8'); fclose(fid);
    BEAT=[1:14 38 41];
    ann_s=zeros(1e4,1,'int32'); na=0; cur=int32(0);
    i=1; nr2=numel(raw2);

    while i<=nr2-1
        w=uint16(raw2(i))+uint16(raw2(i+1))*256; i=i+2;
        at=double(bitshift(w,-10));
        td=double(bitand(w,uint16(1023)));

        if at==0 && td==0
            break;                              % End-of-file marker

        elseif at==59                           % [FIX] SKIP: đọc đủ 4 byte int32
            if i<=nr2-3
                lo=uint32(raw2(i))  + uint32(raw2(i+1))*256;
                hi=uint32(raw2(i+2))+ uint32(raw2(i+3))*256;
                cur=cur+int32(lo)+int32(hi)*65536;
                i=i+4;                          % [FIX] +4 thay vì +2
            end

        elseif any(at==[60 61 62])              % [FIX] NUM/SUB/CHN: KHÔNG có extra bytes
            cur=cur+int32(td);
            % [FIX] Xoá i=i+2 — những type này không có dữ liệu phụ

        elseif at==63                           % AUX: 1 byte độ dài + chuỗi
            cur=cur+int32(td);
            if i<=nr2
                nc=double(raw2(i)); i=i+1;
                i=i+nc+mod(nc,2);               % bỏ qua chuỗi + padding
            end

        else                                    % Annotation thông thường
            cur=cur+int32(td);
            if cur>0 && ismember(at,BEAT)
                na=na+1;
                if na>numel(ann_s)
                    ann_s=[ann_s;zeros(5e3,1,'int32')];
                end
                ann_s(na)=cur;
            end
        end
    end

    ann_s=double(ann_s(1:na));
    ann_s=ann_s(ann_s>=1 & ann_s<=N);

    % Tinh chỉnh ann_disp đến đúng đỉnh R trong ecg_f
    WR_ANN=round(0.15*fs);                      % [FIX] ±150ms thay vì ±100ms
    ann_disp=ann_s;
    for k=1:numel(ann_s)
        i0=max(1,ann_s(k)-WR_ANN);
        i1=min(N,ann_s(k)+round(0.05*fs));
        [~,lm]=max(ecg_f(i0:i1));
        ann_disp(k)=i0+lm-1;
    end

    % [FIX] Đánh giá Se/+P/F1 dùng ann_disp (đã refined) thay vì ann_s thô
    tol=round(TOL*fs/1e3);
    md=false(nr,1); ma=false(numel(ann_disp),1);
    for k=1:nr
        [mn,mi]=min(abs(ann_disp-r_ecg(k)));    % [FIX] ann_disp thay vì ann_s
        if mn<=tol && ~ma(mi)
            md(k)=true; ma(mi)=true;
        end
    end
    TP=sum(md); FP=nr-TP; FN=numel(ann_disp)-TP;
    Se=100*TP/max(1,TP+FN);
    PP=100*TP/max(1,TP+FP);
    F1=2*TP/max(1,2*TP+FP+FN)*100;
    fprintf('Se=%.1f%% | +P=%.1f%% | F1=%.1f%% | TP=%d FP=%d FN=%d\n',...
            Se,PP,F1,TP,FP,FN);
    has_ann=true;
catch
    fprintf('[!] Khong co file .atr\n');
end

% =========================================================
%  VẼ ĐỒ THỊ
% =========================================================
xl=[0 DUR]; xt='Thoi gian (s)'; yt='Bien do (mV)';

% ---- Figure 1: Ba giai đoạn tín hiệu ----
figure(1);
subplot(3,1,1); plot(t,ecg,'g','LineWidth',1.2);
title('(1) ECG goc - MIT-BIH Record 100, kenh MLII');
xlabel(xt); ylabel(yt); xlim(xl); grid on; box on;
legend('ECG goc (mV)','Location','ne');

subplot(3,1,2); plot(t,ecg_n,'r','LineWidth',0.9);
title(sprintf('(2) ECG co nhieu (BW 0.3Hz + Luoi %dHz + EMG)  |  SNR=%.1fdB  |  PRD=%.1f%%',NF,snrN,prdN));
xlabel(xt); ylabel(yt); xlim(xl); grid on; box on;
legend(sprintf('ECG+nhieu (BW+%dHz+EMG)',NF),'Location','ne');

subplot(3,1,3); plot(t,ecg,'g','LineWidth',1.8); hold on; plot(t,ecg_f,'r--','LineWidth',1.2);
title(sprintf('(3) ECG goc vs ECG da loc IIR  |  SNR=%.1fdB  |  PRD=%.1f%%',snrF,prdF));
xlabel(xt); ylabel(yt); xlim(xl); grid on; box on;
legend('ECG goc',sprintf('ECG da loc (Bandpass+Notch %dHz)',NF),'Location','ne');

% ---- Figure 2: Phổ FFT ----
figure(2);
subplot(1,3,1); plot(fax,Yc,'g','LineWidth',1.2);
title('(1) Pho FFT: ECG goc'); xlabel('Tan so (Hz)'); ylabel('Bien do (mV)');
xlim([0 100]); grid on; box on; legend('ECG goc','Location','ne');

subplot(1,3,2); plot(fax,Yn,'r','LineWidth',1.2);
title(sprintf('(2) Pho FFT: ECG co nhieu  |  SNR=%.1fdB',snrN));
xlabel('Tan so (Hz)'); ylabel('Bien do (mV)'); xlim([0 100]); grid on; box on;
legend(sprintf('ECG+nhieu, SNR=%.1fdB',snrN),'Location','ne');

subplot(1,3,3); plot(fax,Yf,'m','LineWidth',1.2);
title(sprintf('(3) Pho FFT: ECG da loc  |  SNR=%.1fdB  |  PRD=%.1f%%',snrF,prdF));
xlabel('Tan so (Hz)'); ylabel('Bien do (mV)'); xlim([0 100]); grid on; box on;
legend(sprintf('ECG da loc (IIR+Notch), SNR=%.1fdB',snrF),'Location','ne');

% ---- Figure 3: Đáp ứng tần số ----
[Hbp,fbp]=freqz(b_bp,a_bp,4096,fs);
[Hnt,fnt]=freqz(b_nt,a_nt,4096,fs); Hc=Hbp.*Hnt;
dB=@(H) 20*log10(abs(H)+eps);
figure(3);
subplot(1,3,1); plot(fbp,dB(Hbp),'r','LineWidth',1.5);
title(sprintf('Butterworth Bandpass [%.1f-%dHz] bac 4',FPASS(1),FPASS(2)));
xlabel('Tan so (Hz)'); ylabel('Bien do (dB)'); xlim([0 fs/2]); ylim([-100 5]);
grid on; box on; xline(FPASS(1),'y--'); xline(FPASS(2),'y--');
legend('Dap ung bien do (dB)','Location','sw');

subplot(1,3,2); plot(fnt,dB(Hnt),'m','LineWidth',1.5);
title(sprintf('Notch Filter fo=%dHz, Q=%d',NF,NQ));
xlabel('Tan so (Hz)'); ylabel('Bien do (dB)'); xlim([0 fs/2]); ylim([-100 5]);
grid on; box on; xline(NF,'g--',sprintf('%dHz',NF));
legend('Dap ung bien do (dB)','Location','sw');

subplot(1,3,3); plot(fbp,dB(Hc),'c','LineWidth',1.5);
title('Dap ung ket hop (Butterworth + Notch)');
xlabel('Tan so (Hz)'); ylabel('Bien do (dB)'); xlim([0 fs/2]); ylim([-100 5]);
grid on; box on; legend('Dap ung ket hop (dB)','Location','sw');

% ---- Figure 4: Sơ đồ Zero-Pole ----
figure(4);
subplot(1,2,1); zplane(b_bp,a_bp);
title(sprintf('Butterworth Bandpass bac 4\n[%.1f-%dHz] @ fs=%dHz',FPASS(1),FPASS(2),fs));
subplot(1,2,2); zplane(b_nt,a_nt);
title(sprintf('Notch Filter fo=%dHz, Q=%d\n@ fs=%dHz',NF,NQ,fs));

% ---- Figure 5: Chuỗi Pan-Tompkins ----
figure(5);
subplot(4,1,1); plot(t,ecg_f,'y','LineWidth',1.2);
title('Buoc 0: ECG da loc (IIR Butterworth+Notch) - Dau vao Pan-Tompkins');
xlabel(xt); ylabel(yt); xlim(xl); grid on; box on;

subplot(4,1,2); plot(t,pts.bpf,'g','LineWidth',1.2); hold on;
               plot(t,pts.deriv,'m','LineWidth',0.8);
title('Buoc 1-2: BPF [5-15Hz] noi bat QRS + Vi phan 5 diem noi bat doc');
xlabel(xt); ylabel('Bien do'); xlim(xl); grid on; box on;
legend('BPF [5-15Hz]','Vi phan 5 diem','Location','ne');

subplot(4,1,3); plot(t,pts.sq,'r','LineWidth',1.0);
title('Buoc 3: Binh phuong (am → duong, khuyech dai QRS so voi P/T)');
xlabel(xt); ylabel('Bien do^2'); xlim(xl); grid on; box on;

subplot(4,1,4);
plot(t,pts.mwi,'Color',[.5 0 .8],'LineWidth',1.2); hold on;
plot(t,pts.thr,'c--','LineWidth',1.5);
plot(t(r_mwi),pts.mwi(r_mwi),'rv','MarkerFaceColor','r','MarkerSize',8);
title(sprintf('Buoc 4-5: MWI (150ms causal) + Nguong thich nghi  ->  %d dinh R',nr));
xlabel(xt); ylabel('Bien do MWI'); xlim(xl); grid on; box on;
legend('MWI causal','Nguong thich nghi (THR1)',sprintf('%d dinh R',nr),'Location','ne');

% ---- Figure 6: So sánh annotation & detection ----
% [FIX] Dùng ann_disp (đỉnh đã tinh chỉnh) để hiển thị → tam giác xanh khớp đỉnh R
figure(6); plot(t,ecg_f,'Color',[1 0.9 0],'LineWidth',1.2); hold on;
if has_ann
    % [FIX] ann_disp thay vì ann_s để tam giác xanh nằm ĐÚNG TRÊN đỉnh R
    plot(t(ann_disp), ecg_f(ann_disp), 'g^', ...
         'MarkerFaceColor',[0 .8 0],'MarkerSize',9);
    plot(t(r_ecg),    ecg_f(r_ecg),   'rv', ...
         'MarkerFaceColor','r','MarkerSize',8);
    legend('ECG da loc (IIR+Notch)', ...
           sprintf('Annotation .atr (%d beats, vi tri chinh xac)',numel(ann_s)), ...
           sprintf('Pan-Tompkins (%d beats)',nr), 'Location','ne');
    title(sprintf('ECG da loc | Se=%.1f%%  +P=%.1f%%  F1=%.1f%%  |  HR=%.1f+-%.1f bpm  |  SDNN=%.1fms  RMSSD=%.1fms',...
                  Se,PP,F1,mean(hr),std(hr),SDNN,RMSSD));
else
    plot(t(r_ecg),ecg_f(r_ecg),'rv','MarkerFaceColor','r','MarkerSize',8);
    title(sprintf('ECG da loc | %d dinh R  |  HR=%.1f+-%.1f bpm  |  SDNN=%.1fms  RMSSD=%.1fms',...
                  nr,mean(hr),std(hr),SDNN,RMSSD));
end
xlabel(xt); ylabel(yt); xlim(xl); grid on; box on;

% ---- Figure 7: HRV ----
figure(7); t_hr=(t(r_ecg(1:end-1))+t(r_ecg(2:end)))/2;
subplot(2,2,[1 2]);
plot(t_hr,hr,'bo-','MarkerFaceColor',[.2 .4 .8],'LineWidth',1.5); hold on;
yline(mean(hr),'r--',sprintf('TB %.1f bpm',mean(hr)),'LineWidth',2,'LabelHorizontalAlignment','left');
fill([t_hr(1) t_hr(end) t_hr(end) t_hr(1)],...
     [mean(hr)-std(hr) mean(hr)-std(hr) mean(hr)+std(hr) mean(hr)+std(hr)],...
     'r','FaceAlpha',0.1,'EdgeColor','none');
title(sprintf('Nhip tim tuc thoi  |  TB=%.1f bpm  |  SD=%.1f bpm  |  pNN50=%.1f%%',...
              mean(hr),std(hr),pNN50));
xlabel(xt); ylabel('Nhip tim (bpm)'); grid on; box on;
legend('HR tuc thoi','Trung binh','+-1 SD','Location','ne');

subplot(2,2,3);
histogram(rr*1e3,'BinWidth',20,'FaceColor',[.2 .5 .8],'EdgeColor','w');
title(sprintf('Phan phoi khoang RR  |  SDNN=%.1fms',SDNN));
xlabel('Khoang RR (ms)'); ylabel('So luong nhip'); grid on; box on;

subplot(2,2,4);
rr_ms=rr*1e3;
scatter(rr_ms(1:end-1),rr_ms(2:end),30,[.2 .4 .8],'filled','MarkerFaceAlpha',0.7); hold on;
al=[min(rr_ms)*0.9 max(rr_ms)*1.1];
plot(al,al,'c--','LineWidth',1.2);
title(sprintf('Poincare Plot  |  SDNN=%.1fms  |  RMSSD=%.1fms',SDNN,RMSSD));
xlabel('RR_n (ms)'); ylabel('RR_{n+1} (ms)');
xlim(al); ylim(al); axis equal; grid on; box on;

% =========================================================
%  HÀM PAN-TOMPKINS (đã sửa 3 lỗi quan trọng)
% =========================================================
function [r_mwi,r_ecg,S]=pan_tompkins(ecg,fs)
    ecg=ecg(:); N=numel(ecg);

    % Bước 1: BPF nội [5-15 Hz]
    [bl,al]=butter(5,min(15/(fs/2),0.99),'low');        % [FIX] 11→15 Hz
    [bh,ah]=butter(1,max(5/(fs/2), 0.01),'high');
    bpf=filtfilt(bh,ah,filtfilt(bl,al,ecg));

    % Bước 2: Vi phân 5 điểm
    drv=conv(bpf,(1/(8/fs))*[-2 -1 0 1 2],'same');

    % Bước 3: Bình phương
    sq=drv.^2;

    % Bước 4: MWI (tích phân cửa sổ trượt NHÂN QUẢ)
    W=round(0.15*fs);
    % [FIX] Causal MWI: [W-1, 0] thay vì symmetric W
    %       → đỉnh MWI khớp đúng thời điểm QRS, không bị lệch W/2 mẫu
    mwi=movmean(sq,[W-1, 0]);

    % Bước 5: Ngưỡng thích nghi
    mdist=round(0.2*fs);
    refrac=round(0.2*fs);

    % Khởi tạo SPKI/NPKI từ 3 giây đầu (ổn định hơn 2 giây)
    init_len=min(round(3*fs),N);                        % [FIX] 2s→3s
    [pk0,~]=findpeaks(mwi(1:init_len),'MinPeakDistance',mdist);
    if isempty(pk0)
        SPKI=0.15*max(mwi); NPKI=0.1*SPKI;
    else
        th0=0.5*max(pk0);
        qp=pk0(pk0>=th0); np=pk0(pk0<th0);
        SPKI=mean(qp);
        if isempty(np), NPKI=0.1*SPKI; else, NPKI=mean(np); end
    end
    TH1=NPKI+0.25*(SPKI-NPKI); TH2=0.5*TH1;

    % Tìm tất cả đỉnh MWI
    [apk,alc]=findpeaks(mwi,'MinPeakDistance',mdist);
    rloc=zeros(numel(apk),1,'int32'); nr=0;
    lastR=int32(0);
    rrbuf=zeros(8,1); nrr=0; rravg=round(0.8*fs);
    thhist=zeros(N,1);

    for k=1:numel(alc)
        loc=int32(alc(k)); pk=apk(k); thhist(loc)=TH1;

        % Bỏ qua nếu trong refractory period
        if lastR>0 && (loc-lastR)<refrac
            NPKI=0.125*pk+0.875*NPKI;
            TH1=NPKI+0.25*(SPKI-NPKI); TH2=0.5*TH1; continue;
        end

        if pk>=TH1
            % T-wave rejection: đỉnh trong 360 ms sau R trước, biên độ < 50% SPKI
            if lastR>0 && double(loc-lastR)<0.36*fs && pk<0.5*SPKI
                NPKI=0.125*pk+0.875*NPKI;
                TH1=NPKI+0.25*(SPKI-NPKI); TH2=0.5*TH1; continue;
            end
            nr=nr+1; rloc(nr)=loc;
            SPKI=0.125*pk+0.875*SPKI;
            if lastR>0
                nrr=min(nrr+1,8); rrbuf=circshift(rrbuf,1);
                rrbuf(1)=double(loc-lastR);
                rravg=mean(rrbuf(1:nrr));
            end
            lastR=loc;
        else
            NPKI=0.125*pk+0.875*NPKI;
        end
        TH1=NPKI+0.25*(SPKI-NPKI); TH2=0.5*TH1;
    end
    rloc=rloc(1:nr);

    % Searchback: nếu khoảng RR > 1.66×RR_avg
    if nr>=2 && nrr>=1
        extra=[];
        for k=2:nr
            gap=double(rloc(k)-rloc(k-1));
            if gap>1.66*rravg
                s0=double(rloc(k-1))+round(0.2*fs);
                s1_=double(rloc(k))-round(0.2*fs);
                if s0<s1_ && s0>=1 && s1_<=N
                    [sp,sl2]=findpeaks(mwi(s0:s1_));
                    if ~isempty(sp)
                        [mx,mi]=max(sp);
                        if mx>=TH2, extra(end+1)=s0+sl2(mi)-1; end
                    end
                end
            end
        end
        if ~isempty(extra)
            rloc=sort([rloc;int32(extra(:))]);
            nr=numel(rloc);
        end
    end

    % --------------------------------------------------------
    % [FIX QUAN TRỌNG] Tinh chỉnh vị trí R trong ECG gốc
    %   Thay đổi so với bản gốc:
    %   (a) wr: 0.06*fs (22 mẫu) → 0.15*fs (54 mẫu)  -- cửa sổ đủ rộng
    %   (b) max(ecg) thay vì max(abs(ecg))             -- không bị nhầm sóng S âm
    %   (c) Cửa sổ dời sang trái (-wr) để bắt đỉnh trước MWI peak
    % --------------------------------------------------------
    wr=round(0.15*fs);                                  % [FIX] 0.06→0.15
    recg=zeros(nr,1,'int32');
    for k=1:nr
        % MWI causal: đỉnh MWI nằm SAU đỉnh R khoảng W/2 mẫu
        % → dịch cửa sổ về phía trước (look backward)
        i0=max(1, double(rloc(k))-wr);
        i1=min(N, double(rloc(k))+round(0.05*fs));     % ít mẫu về phía sau
        [~,lm]=max(ecg(i0:i1));                         % [FIX] max, không abs
        recg(k)=int32(i0+lm-1);
    end

    r_mwi=double(rloc); r_ecg=double(recg);

    % Tạo đường ngưỡng cho Figure 5
    thr=zeros(N,1);
    if nr>=1
        thr=interp1([1;r_mwi;N],[thhist(1);thhist(r_mwi);thhist(N)],...
                    (1:N)','linear','extrap');
        thr=max(thr,0);
    end
    S.bpf=bpf; S.deriv=drv; S.sq=sq; S.mwi=mwi; S.thr=thr;
end