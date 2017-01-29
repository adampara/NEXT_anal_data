
import matplotlib.pyplot as plt
from ROOT import gRandom, TCanvas, TH1F, TFile, TTree,TH2F,gDirectory, TF1,TF2,gStyle
from scipy import interpolate
import ROOT
import sys
from math import exp
from plot_histogram import plot_histogram
from smooth_wave import smooth_wave
#import peakutils
from peakutils.plot import plot as pplot
# from peakdetect import peakdetect
from detect_peaks import detect_peaks
from Browse_ROOT import Browse_ROOT
import pickle

def most_common(lst):
    return max(set(lst), key=lst.count)

def ac_coupl(wf,tau):
    att = exp(-tau)

    bs = [0 for i in range(len(wf))]
    wf_corr = [0 for i in range(len(wf))]

    bss = 0
    for i in range(1,len(wf)):
##         wf_corr[i] = wf[i]-bs[i-1]*att
##         bs[i] = bs[i-1]*att - wf_corr[i]*tau
        wf_corr[i] = wf[i]-bss*att
        bss = bss*att - wf_corr[i]*tau

    return wf_corr,bs

def GetKeyNames ( self) :
    return [name.GetName() for name in self.GetListofKeys()]

TFile.GetKeyNames = GetKeyNames

def ldir(dirname,MyFile, Sel_Dir='',debug=False):
    """
    return list af all histgrams in a root file
    """
    
    hh= []      # list of histograms

    savdir = gDirectory.GetPathStatic()
    savd = gDirectory
    if debug:
        print '   initial directory  ', savd, '        ', savdir
    
    lk1 = dirname.GetListOfKeys()    

    for key in lk1:

        if debug:
            print '  =======  Key found  ',key.GetName(),key.GetClassName()
        if(key.GetClassName() == 'TH1F'):
            #  check if this directory is selected
            if len(Sel_Dir) < 1 or (
                len(Sel_Dir)>1 and savdir.split('/')[-1] == Sel_Dir):
                print savdir, Sel_Dir
                hh.append(key)

                # below is an example of accessing a tree.
        if  key.GetClassName() == 'TTree':
             print '    Class = ', key.GetClassName()
             tree = MyFile.Get(key.GetName())
             print tree

             branches =  tree.GetListOfBranches()
             print '  branches ',branches
             print branches.GetEntries()

             for i in range(0,branches.GetEntries()):
                 branch = branches.At(i)
                 print branch
                 print branch.GetName()

             leaves = tree.GetListOfLeaves()
             print '  leaves list of entries'
             print leaves
             print leaves.GetEntries()

             for i in range(0,leaves.GetEntries()):
                 leaf = leaves.At(i)
                 print leaf
                 print leaf.GetName()
                 #    GetLen is the number of elements
                 #    it is 0 if it is specified in another leaf
                 print leaf.GetLen()

             tree.GetEntry(15)
             print ' Number of waveform points '
             print leaves.At(0).GetValue()
             np = int(leaves.At(0).GetValue())
             print np
             for i in range (0,np):
                 print i, leaves.At(1).GetValue(i)


        if  key.GetClassName() == 'TDirectoryFile':

            dirname.cd(key.GetName())
            
            if debug:
                print 'This key is a directory,', key.GetName()
                print ' Current directory is '
            gDirectory.pwd()

            dirn = gDirectory
            hh.extend(ldir(dirn,MyFile,Sel_Dir=Sel_Dir))

            gDirectory.cd('..')
            gDirectory.pwd()

            savd.cd()
            gDirectory.pwd()
            if debug:
                print savdir
            gDirectory.cd(savdir)
            gDirectory.pwd()

            continue

            if debug:
                print '   after the loop, savdir ' ,savdir
                print '   after the loop, savd ' ,savd
                print '     end of a listing of a directory  '

    return hh

import numpy as np
def parse(lines):

    mon = {'Jan': '01', 'Feb': '02', 'Mar': '03',
           'Apr': '04', 'May': '05', 'Jun': '06',
           'Jul': '07', 'Aug': '08', 'Sep': '09',
           'Oct': '10', 'Nov': '11', 'Dec': '12'}
    wform = []
    acq_rec = []
    acq_no = -1
    tim = []
    ampl = []

    for i in range(len(lines)):
        #print 'line number ',i
        if i<10:
            print i,'  ',lines[i]
        if i%11 == 0:

            wf = lines[i].split()[1]

        if i%11 == 1:
            acq = lines[i].split()[1]

        if i%11 == 2:
            dt = float(lines[i].split()[1])

        if i%11 == 9:
            if int(wf) < 10:
                print (lines[i])
            tok = lines[i][1:len(lines[i])-3].split()
            tso = lines[i][1:len(lines[i])-3]
            tsa = tok[2]+'-'+tok[0]+'-'+mon[tok[1]]+'T'+tok[3]+tok[4]+tok[5]+'Z'
            tsb = np.datetime64(tsa)
            if i == 9:
                T0 = tsb
            tsc = tsb - T0

            ts = int(str(tsc).split()[0])*1.e-9

        if i%11 == 10:
            entry = (wf,acq,ts)


            if acq_no != acq:
                #  new acquisition, if not the first one, store the previous one
                print 'old acquisition ',acq_no, '  new acquisition ',acq
                if acq_no == -1:
                    print 'the first acquisition'
                    acq_rec.append(entry)
                else:
                    print 'the next acquisition'
                    wform.append(acq_rec)
                    acq_rec = []
                    acq_rec.append(entry)
                acq_no = acq
            else:
                acq_rec.append(entry)
            if i == len(lines)-1 :
                print 'last frame'
                wform.append(acq_rec)

            #  frame

            aa = lines[i].split()
            for i in range(len(aa)):
                tt = ts + i*dt
                tim.append(tt)
                ampl.append(float(aa[i]))

    return wform,tim,ampl

    
def process_SiPM(h, debug=False):
    """
    process a single SiPM waveform from a root histograam
    """
    
    if debug:
        print 'process histogram ', h.GetName()
    
    # extract information about the histogram
    nbin = h.GetNbinsX()
    xc = np.zeros(nbin)
    cont = np.zeros(nbin)
    for bin in range(nbin):
        xc[bin] = h.GetBinCenter(bin+1)
        cont[bin] = h.GetBinContent(bin+1)
    delta = xc[1] - xc[0]       # time binning

    offset = 1      # to cut the backgorud
    baseline = most_common(cont.tolist()) 
    if debug:
        print 'overall baseline ', baseline

    cont_bas = cont  -  baseline  # subtract baseline, invert the waveform
    cont_b = np.zeros(nbin)
    for i in range(len(cont_bas)):
        cont_b[i] = max(0,cont_bas[i]-offset)
        
    tfront = int(590/delta)
    t_end_sig = int(630/delta)
    t_end = int(1200/delta)
    
    q_front = np.sum(cont_b[:tfront])
    q_post = np.sum(cont_b[t_end_sig:t_end])
    q_sig = np.sum(cont_b[tfront:t_end_sig])
    
    if debug:
        print '   qfront, qtail ',q_front, q_post
        plt.plot(xc,cont_b)
        plt.show()
    return q_front,q_post,q_sig
    

    
def process_histogram(h, debug=False):
    """
    process a single waveform froma root histograam
    """
    
    if debug:
        print 'process histogram ', h.GetName()
    
    # extract information about the histogram
    nbin = h.GetNbinsX()
    xc = np.zeros(nbin)
    cont = np.zeros(nbin)
    for bin in range(nbin):
        xc[bin] = h.GetBinCenter(bin+1)
        cont[bin] = h.GetBinContent(bin+1)
    delta = xc[1] - xc[0]       # time binning

    baseline = most_common(cont.tolist())
    if debug:
        print 'overall baseline ', baseline

    cont_b = baseline - cont    # subtract baseline, invert the waveform
    nsmooth = 5
    cont_smooth = smooth_wave(cont_b,nsmooth)  # smooth out the waveform

    # select the waveforms
    tfront = int(590/delta)
    t_end_sig = int(650/delta)
    t_end = int(1200/delta)
    t_end_sig_short = int(620/delta)
    good_front_f = []
    good_post_f = []
    good_front_f_Q = []
    good_post_f_Q = []
    
    s1_thr = 20.
    mpd_s1 = 50        # minimum peak distance
    ind_front = detect_peaks(cont_smooth[:tfront], mph = s1_thr, mpd = mpd_s1)
    s2_thr = 200.
    mpd_s2 = 200
    ind_sig = detect_peaks(cont_smooth[tfront:t_end_sig], mph = s2_thr, 
                           mpd = mpd_s2) + tfront
    end_thr = 20.
    mpd_end = 200
    ind_end = detect_peaks(cont_smooth[t_end_sig:t_end], mph = end_thr, 
                           mpd = mpd_end) + t_end_sig
    sig = np.sum(cont_smooth[tfront:t_end_sig]) 
    sig_sh = np.sum(cont_smooth[tfront:t_end_sig_short])
    
    clean = False
    
    #   select short S2 peak for ad hoc analysis
    
    if len(ind_front) == 1 and len(ind_sig)>0 and len(ind_end) == 0:
        clean = True
        

    if clean:
        pe_thr = 2.0
        pe_width = 0.2
        pe_len = int(0.5 * pe_width/delta)
        
        spe_front = detect_peaks(cont_smooth[:tfront], mph = pe_thr, 
                                 mpd = pe_len)
        spe_post  = detect_peaks(cont_smooth[t_end_sig:t_end], mph = pe_thr, 
                                 mpd = pe_len)
    
        (good_post, bad_post) = select_peaks(cont_smooth[t_end_sig:t_end],
                    pe_thr,pe_len,spe_post, debug = False)
        (good_front, bad_front) = select_peaks(cont_smooth[:tfront],
                    pe_thr,pe_len,spe_front, debug = False) 
        good_post_f = [ x + t_end_sig for x in good_post]
        good_front_f = list(set(good_front)-set(ind_front))  # subtract S1
        
        for i in good_front_f:
            good_front_f_Q.append(np.sum(cont_smooth[i-pe_len:i+pe_len]))
        for i in good_post_f:
            good_post_f_Q.append(np.sum(cont_smooth[i-pe_len:i+pe_len]))
        

        
    if debug and clean:
        print ' Q signal ',sig
        print ' peaks found, front', xc[ind_front]
        print ' pe peaks, front ', xc[good_front_f]
        print ' pe Q, front ', good_front_f_Q
        print ' peaks found, signal', xc[ind_sig]
        print ' peaks found, end', xc[ind_end]
        print ' pe Q, end ', good_post_f_Q
        
        #plt.plot(xc,cont_smooth)
        pplot(xc,cont_smooth, ind_front)
        pplot(xc,cont_smooth, ind_sig)
        pplot(xc,cont_smooth, ind_end)
        plt.show()

        # show accepted peaks
        pplot(xc,cont_smooth,good_front_f+good_post_f)
        plt.show()
    
    return (sig, xc[good_front_f], xc[good_post_f],
            cont_smooth[good_front_f], cont_smooth[good_post_f],
            good_front_f_Q, good_post_f_Q, clean)
  
def select_peaks(wf, thresh, wpeak, lpeak, debug = False):
    """
    find a list of peaks with given width
    """
    good_peaks = []
    bad_peaks = []
    if debug:
        print 'expected peak width ', wpeak 
        
    for peak in lpeak:

        acc = accept_peak(wf,peak,thresh,wpeak, debug = debug)

        if acc:
            good_peaks.append(peak)
        else:
            bad_peaks.append(peak)
    
    return good_peaks, bad_peaks
    
def accept_peak(wf,peak,thresh,wpeak, debug = False):
    """
    accpet a single pe peak
    """

    thrval = 0.7
    window = 10     #  window to get local baseline
    
    acc = False
    if debug:
        print peak
    if peak > wpeak + window:
        if debug:
            print ' beginning ok ', peak
        if len(wf)-peak - window > wpeak:
            if debug:
                print 'end ok ',peak
                print 'height, left', wf[peak] - np.mean(wf[peak-wpeak-window:peak-wpeak])
            if wf[peak] - np.mean(wf[peak-wpeak-window:peak-wpeak]) > thrval*thresh:
                if debug:
                    print ' left side OK', peak
                    print ' height, right', wf[peak] - np.mean(wf[peak+wpeak:peak+wpeak+window])
                if  wf[peak] - np.mean(wf[peak+wpeak:peak+wpeak+window]) > thrval*thresh:
                    if debug:
                        print 'right side ok ', peak
                    if np.mean(wf[peak-2:peak+2]) > 1.0*thresh:     
                        if debug:
                            print ' waveform at the peak', wf[peak-2:peak+2]
                            print '  mean around the peak', np.mean(wf[peak-2:peak+2])
                        acc = True
    if debug:
        print 'accept = ', acc
        ibeg = max(0, peak - 100)
        iend = min(len(wf),peak + 100)
        plt.plot(wf[ibeg:iend])
        plt.plot(peak-ibeg,wf[peak],'*')
        plt.show()    
    return acc
    
    """Detect peaks in data based on their amplitude and other features."""

#from __future__ import division, print_function

import os
import numpy as np
import glob

def Book_Hist_SiPM():
    """
    Book histograms for single pe analysis
    """
    
    lhist = []
    
    h1 = TH1F('q_SiPM_f', 'SiPM charge front', 100, 1, 0)
    lhist.append(h1)
    h2 = TH1F('q_SiPM_p', 'SiPM charge post', 100, 1, 0)
    lhist.append(h2)
    h3 = TH1F('q_SiPM_diff', 'SiPM charge post - front', 100, 1, 0)
    lhist.append(h3)
    h4 = TH1F('q_SiPM_sig', 'SiPM signal', 100, 1, 0)
    lhist.append(h4)
    
    return (lhist,h1,h2, h3, h4)
    
filedir = 'data/data_Vb_26_5_Vth_5_380nm_50nW_5_Ch9-2016-12-13_19-49-25/'
filenam = 'data_Vb_26_5_Vth_5_380nm_50nW_5_Ch9-2016-12-13_19-49-25.dat'
filedir = 'data/noise_Vb_27_5_Vth_10_Ch9-2016-12-16_16-39-37/'
filedir = '/Volumes/Seagate Backup Plus Drive/DATA/Hamamatsu/S13360-3050CS/AC/from_45C_to_m175/Miteq20MHz/DCR/03_04_2016/45/'
filedir = 'data/DCR/03_04_2016/-135/'
filedir = 'data/DCR/03_04_2016/'
# filelist = os.listdir(filedir)

an_fil = 'run_2983'
an_fil = 'run_3112'
filelist = glob.glob(an_fil+'*ALL*root*')
print filelist

debug = False
#file0 = "run_3112_000_0_100.root"
# file0 = "run_2983_000_0_100.root"
#(lhist,h1,h2,h3,h4,h5,h6,h7,h8,h9) = Book_Hist_Pe()
#filelist = [file0]

cle_ev = ['216', '217', '768', '1146', '213', '1302', '1143', '1262', '1308',
            '1269', '218', '957', '1147', '24', '20', '1076', '1072', '29', 
            '94', '542', '406', '1255', '403', '1468', '401', '545', '933', 
            '812', '930', '571', '408', '1097', '1158', '719', '1691', '1153', 
            '1309', '421', '1492', '260', '267', '1265', '265', '1497', '1702', 
            '1625', '1652', '892', '417', '1258', '372', '776', '1086', '988', 
            '1321', '1082', '597', '596', '599', '194', '629', '196', '1124', 
            '190', '985', '270', '88', '89', '397', '396', '82', '83', '1080', 
            '87', '84', '797', '773', '790', '799', '1133', '1254', '725', 
            '918', '368', '1251', '366', '364', '913', '361', '1339', '443', 
            '1131', '1130', '1137', '446', '1135', '244', '246', '1514', '100', 
            '1298', '1436', '785', '787', '781', '904', '1642', '896', '1658', 
            '30', '789', '1438', '34', '1521', '926', '1247', '623', '1683', 
            '1323', '573', '451', '452', '1343', '455', '456', '574', '60', 
            '1346', '259', '179', '256', '1504', '745', '976', '975', '931', 
            '1634', '1081', '1427', '1637', '1429', '747', '1301', '731', 
            '748', '180', '184', '6', '1469', '980', '568', '1123', '634', 
            '544', '560', '561', '400', '98', '90', '222', '96', '1623', 
            '1483', '1466', '1627', '1626', '17', '1629', '1628', '727', 
            '273', '1452', '1453', '1487', '1462', '934', '729', '1165', '607', 
            '558', '746', '602', '555', '554', '557', '556', '1640', '1161', 
            '1439', '231', '756', '1681', '1325', '46', '950', '953', '1613', 
            '40', '41', '1696', '1110', '1113', '1299', '1505', '1690', '1440', 
            '9', '1293', '1296', '1449', '891', '203', '619', '207', '206', 
            '1522', '611', '1274', '1276', '75', '923', '72', '1636', '79', 
            '792', '1285', '1284', '1286', '1281', '1312', '942', '1607', 
            '1472', '1616', '809', '1684', '1689', '572', '358']
clean_ev = [int(x) for x in cle_ev]
print clean_ev

gDirectory.pwd()
dirhome = '/'
print 'home directory', dirhome
gDirectory.cd(dirhome)

(lhist,h1,h2, h3, h4) = Book_Hist_SiPM()
pfile = open('sipm.pkl', 'rb')
SiPM_ID = pickle.load(pfile)
print SiPM_ID

sipm_diff = {}
for sipm in SiPM_ID:
    handle = 'sipm_' + sipm
    title = 'Tail - front, sipm ' + sipm
    gDirectory.cd(dirhome)
    sipm_diff[sipm] =  TH1F(handle, title, 100, 0, 2000)
    lhist.append(sipm_diff[sipm])
    
ifile = 0
for finp in filelist:
    ifile += 1
    if ifile > 20:
        continue
    print 'ifile ',ifile
    MyFile = TFile(finp)
    
    gDirectory.pwd()
    dir = gDirectory
    if debug:    
        print 'file directory  ',dir
     
    if debug:   
        # print hists
        print 'list of histograms'
    
    event_hist = []
    last_ev = '0'
    n_anal = 0
    max_anal = 9990000000
    
    histSi = ldir(dir,MyFile,Sel_Dir='histSi')

    #histP = ldir(dir,MyFile,Sel_Dir='histP')
    histD = []
    #histD = ldir(dir,MyFile,Sel_Dir='histD')   
    last_ev = '-1'
    for h in histSi:

        if n_anal < max_anal:
        
            hh = h.ReadObj()
            tok = h.GetName().split('_')
            SiPM_ID = int(tok[1])
            evno = int(tok[4])
            #print 'ev number ', evno
            if evno in clean_ev:
                (q_front,q_tail, q_sig) = process_SiPM(hh, debug=False)
                h1.Fill(q_front)
                h2.Fill(q_tail)
                h3.Fill(q_tail - q_front)
                h4.Fill(q_sig)
                sipm_diff[tok[1]].Fill(q_tail)
                n_anal +=1
                print ' clean event ', evno
                if debug:
                     print '  waveform ', hh.GetName()
                     print 'event number ',evno, 'SiPM', SiPM_ID 
                     print '  front, tail' , q_front, q_tail
                     plot_histogram(hh)
            
                if tok[4] != last_ev:
                    #   a new event  
                    if last_ev != '-1':
                        event_hist.append(ev_hist)
                    last_ev = tok[4]
                    ev_hist = hh.Clone()
                    ev_hist.SetName('event_'+tok[4])
            
                else:
                    ev_hist.Add(hh)
        
                 
    for h in histD:
        n_anal += 1
        if n_anal < max_anal:
            hh = h.ReadObj()
            (sig, time_front, time_post, ampl_front, ampl_post,
                    Q_front, Q_post, clean) = process_histogram(hh, debug=False)
            
            if clean:
                h1.Fill(len(time_front))      
                for tf,af,qf in zip(time_front, ampl_front, Q_front):       
                    h3.Fill(qf)
                    h5.Fill(tf)
                    h7.Fill(af)
                h2.Fill(len(time_post))      
                for tp,ap,qp in zip(time_post, ampl_post, Q_post): 
                    h4.Fill(qp)
                    h6.Fill(tp)
                    h8.Fill(ap)
                h9.Fill(sig)
            
            #    -----    event -based analysis    --------
            #   sum up all waveforms of a given event
            
            tok = h.GetName().split('_')
            if tok[4] != last_ev:
                #   a new event  
                if last_ev != '0':
                    event_hist.append(ev_hist)
                last_ev = tok[4]
                ev_hist = hh.Clone()
                ev_hist.SetName('event_'+tok[4])
        
            else:
                ev_hist.Add(hh)
    
    #  addpend the last event to the list
    # event_hist.append(ev_hist)

hist_file = an_fil+'SiPM.hst'
# hist_file = 'run2983.hst'

fr = TFile(hist_file, 'recreate')
for h in lhist:
    h.Write()
fr.Close()

Browse_ROOT()

exit()  
for h in event_hist:
    print h.GetName()
    plot_histogram(h)

