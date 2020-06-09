# coding: utf-8

# Import computational modules
import json
import os
from sys import argv
import time
import pprint

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.signal import convolve2d
from skimage.draw import circle
from scipy.interpolate import RectBivariateSpline

# Surpress any warnings
import warnings
warnings.filterwarnings(action='ignore')


def main(filename, interpol_step=None, particle_no=None):
    iterator(filename, None, second_pass, interpol_step, particle_no)


class Log(dict):
    '''
    Class that contains the info from the log file, irrespective of its version.
    To load log file use Log.load($LOG_FILE_PATH$). Class will make dict out of
    it. To write the instance to file as a new log file run:

        self.write_log($NEW_LOG_FILE_PATH$, type='new')

    To implement:
        - set the dict values as immutable properties
        - make the dict values accessbile with self.mode instead of
          self['detector']['mode']
        - make write_old method
        - make __init__ so it accepts both $LOG_FILE_PATH$ and **kwargs.
          Currently only **kwargs are supported
        - tests and Error raises for weird values (e.g. mode = 2.4)
    '''
    @classmethod
    def load_log(cls, log_file_path):
        # TODO maybe dict instead of if
        is_new = cls.is_new(log_file_path)

        if is_new:
            return cls.load_new_log(log_file_path)

        elif not is_new:
            return cls.load_old_log(log_file_path)

    @classmethod
    def load_new_log(cls, log_file_path):
        with open(log_file_path, 'r') as log_file:
            log_dict = json.load(log_file)
            # pprint(log_dict)
        return cls(**log_dict)

    @classmethod
    def load_old_log(cls, log_file_path):
        log_dict = {'detector': {}}

        with open(log_file_path, 'r') as log_file:
            detector_num = int(log_file_path.strip('.log').split('-D')[-1])

            log_dict['filename'] = os.path.abspath(log_file_path)

            for i, line in enumerate(log_file):
                if i == 0:
                    log_dict['timestamp'] = line.strip('\n')
                elif i == 1:
                    continue
                elif i == 2:

                    mode = int(line.split(' ')[-1])
                    log_dict['detector']['mode'] = mode
                elif i == 3:
                    width, height = [int(item) for item in line.split('x')]
                    log_dict['detector']['image_size'] = {'width': width,
                                                          'height': height
                                                          }
                elif i == 4:
                    left, top, right, bottom = [
                        int(item) for item in line.replace(',', '').split(' ')[-4:]]
                    log_dict['detector']['ROI'] = {'left': left,
                                                   'top': top,
                                                   'right': right,
                                                   'bottom': bottom
                                                   }
                elif i == 5:
                    fps = float(line.split(' ')[-1])
                    log_dict['detector']['fps'] = fps

                    frames = []
                    for j, line in enumerate(log_file):
                        frames.append([j, int(line), detector_num])
                    log_dict['detector']['frames'] = frames

        return cls(**log_dict)

    @classmethod
    def is_new(cls, log_file_path):
        with open(log_file_path, 'r') as log_file:
            try:
                json.load(log_file)
                return True
            except json.JSONDecodeError:
                return False

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def write_log(self, log_file_path, type=None):
        write_method_dict = {'new': self.write_new_log,
                             'old': self.write_old_log
                             }

        return write_method_dict[type](log_file_path)

    def write_new_log(self, new_log_file_path):
        if os.path.exists(new_log_file_path) and not Log.is_new(new_log_file_path):
            os.rename(new_log_file_path, '{}.copy'.format(new_log_file_path))

        with open(new_log_file_path, 'w') as new_log_file:
            output1 = pprint.pformat(self)
            output2 = output1.replace('\'', '"')

            new_log_file.write(output2)

    def write_old_log(self, log_file_path):
        raise NotImplementedError


def load_images(file, first_image=None, last_image=None):
    """
    Function takes file name without the extension, ordinal number of the
    first image, and ordinal number of the last image as the paramaters.
    It returns the array of images between first and last image including the
    first and the last image.
    Note: ordinal number of image is 1 for the first image - just like in Fiji
          Returned array of images has different convention - first image has
          ordinal number of 0.
    """

    log_file = '../../' + file + '.log'
    seq_file = '../../' + file + '.seq'

    if os.path.exists(log_file):
        log_file_dict = Log.load_log(log_file)
        nx = log_file_dict['detector']['image_size']['width']
        ny = log_file_dict['detector']['image_size']['height']
    else:
        nx = 768
        ny = 960

    seq_size = os.path.getsize(seq_file)
    total_images = seq_size // (nx * ny * 2)

    if first_image is None:
        first_image = 0
    else:
        first_image = first_image - 1

    if last_image is None:
        last_image = total_images - 1
    else:
        last_image = last_image - 1

    n_images = last_image - first_image + 1

    with open(seq_file, 'br') as sequence:
        dtype = 'H'
        mode = 'r'
        offset = ny * nx * first_image * 2
        shape = (n_images, ny, nx)
        # Use memmap when RAM is low, and fromfile + seek + reshape when you
        # can load whole .seq file in RAM
        images = np.memmap(sequence,
                           dtype=dtype,
                           mode=mode,
                           offset=offset,
                           shape=shape
                           )

    return images, ny, nx


def load_images_average(file):
    file = os.path.basename(file)
    print(file)
    print('Trying to load {}'.format('averages{}.npy'.format(file)))
    try:
        averages = np.load('../../averages{}.npy'.format(file))
        print('Average image found and loaded!')

    except FileNotFoundError:
        print('Calculating average image...')
        images, *_ = load_images(file)
        averages = np.average(images, axis=0)
        np.save('./Averages/averages{}.npy'.format(file), averages)

    return averages

# For computational feasibility we will look for intruder in a small window (patch) that surrounds it
# Patch is drawn around the center found in previous T.step, thus the size of patch should be big enough
# so that intruder doesn't escape the patch between two subseqeuent images


def average_image_of_patch(T):
    # Index of upper boundary in array
    upper_limit = T.ny - 1
    # Index of right boundary in array
    right_limit = T.nx - 1
    # Patch size is large enough so that intruder doesn't leave the patch in
    # two subsequent frames
    delta_patch = 10
    # delta patch is very important variable - it has to be larger than distance
    # that intruder covers in two subsequent images

    patch_size = int(2 * T.r_real + 2 * delta_patch)

    # Check if patch is outside the image's upper/lower boundary
    if int(T.cy_real) - patch_size // 2 < 0:
        bottom_edge = 0
    # Variable upper limit is passed from log file and it is (vertical resolution - 1)
    elif int(T.cy_real) + patch_size // 2 > upper_limit:
        bottom_edge = upper_limit - patch_size + 1
    else:
        bottom_edge = int(T.cy_real) - patch_size // 2

    # Check if patch is outside the image's left/right boundary
    if int(T.cx_real) - patch_size // 2 < 0:
        left_edge = 0
    # Variable right limit is passed from log file and it is (horizontal resolution - 1)
    elif int(T.cx_real) + patch_size // 2 > right_limit:
        left_edge = right_limit - patch_size + 1
    else:
        left_edge = int(T.cx_real) - patch_size // 2
    # Extract the part of original image that is within the patch
    patch_image = T.image[bottom_edge:bottom_edge + patch_size,
                          left_edge:left_edge + patch_size]
    # Extract the part of averaged image that is within the patch
    patch_averages = T.averages[bottom_edge:bottom_edge + patch_size,
                                left_edge:left_edge + patch_size]
    # Divide patch image with patch average
    averaged_patch = np.divide(patch_image, patch_averages)
    # Get rid of weird values in averaged patch
    averaged_patch[np.isnan(averaged_patch)] = 0

    return averaged_patch, bottom_edge, left_edge


def interpolate(patch, interpolation_step):
    dx1, dy1 = 1, 1
    x1 = np.arange(0, patch.shape[1], dx1)
    y1 = np.arange(0, patch.shape[0], dy1)

    image_interpol = RectBivariateSpline(x=y1,
                                         y=x1,
                                         z=patch,
                                         kx=1,
                                         ky=1
                                         )

    dx2, dy2 = 1 / interpolation_step, 1 / interpolation_step
    x2 = np.arange(0, patch.shape[1], dx2)
    y2 = np.arange(0, patch.shape[0], dy2)

    patch_2 = image_interpol(y2, x2)

    return patch_2


def center_of_darkest_circle(patch, r, interpol=False, cy=None, cx=None, step=None, bottom=None, left=None):
    # Convert r to native Python float type, we need it for round()
    r = float(r)  # .item()
    if interpol == False:
        # This part is run on not interpolated patch
        mask = np.zeros((round(2 * r), round(2 * r)))
        rr, cc = circle((mask.shape[1] - 1) / 2, (mask.shape[0] - 1) / 2, r)
        mask[rr, cc] = 1
        # In case I want weighted function I can use line below instead of 1
        # mask[rr, cc] = 1 - np.sqrt(((cc - (r - 0.5)) ** 2 +
        # (rr - (r - 0.5)) ** 2) / (1.2 * r)**2)

        conv = convolve2d(patch, mask, mode='valid')

        b = np.where(conv == conv.min())
        # Center of circle in the patch is given by:
        cy = (mask.shape[0] - 1) / 2 + b[0][0]
        cx = (mask.shape[1] - 1) / 2 + b[1][0]
    else:
        # This part is run on interpolated patch
        mask = np.zeros((int(2 * r), int(2 * r)))
        rr, cc = circle(r - 0.5, r - 0.5, r)
        mask[rr, cc] = 1
        # New patch center is:
        cx_cen = (cx - left) * step + (step - 1) / 2.
        cy_cen = (cy - bottom) * step + (step - 1) / 2.
        # New patch boundaries are int(c_cen - (r+step/2.)) and int(c_cen + 1 + (r+step/2.))

        if int(cy_cen - (r + step / 2.)) < 0:
            be = 0
            ue = int(1 + 2 * (r + step / 2.))
        else:
            if int(cy_cen + 1 + (r + step / 2.)) > patch.shape[0]:
                ue = patch.shape[0]
                be = ue - int(1 + 2 * (r + step / 2.))
            else:
                be = int(cy_cen - (r + step / 2.))
                ue = int(cy_cen + 1 + (r + step / 2.))

        if int(cx_cen - (r + step / 2.)) < 0:
            le = 0
            re = int(1 + 2 * (r + step / 2.))
        else:
            if int(cx_cen + 1 + (r + step / 2.)) > patch.shape[1]:
                re = patch.shape[1]
                le = re - int(1 + 2 * (r + step / 2.))
            else:
                le = int(cx_cen - (r + step / 2.))
                re = int(cx_cen + 1 + (r + step / 2.))

        patch = patch[be:ue, le:re]

        conv = convolve2d(patch, mask, mode='valid')

        b = np.where(conv == conv.min())
        cy = be + mask.shape[0] / 2 + b[0][0] - 0.5
        cx = le + mask.shape[0] / 2 + b[1][0] - 0.5

    return cy, cx


def open_data_frame(file):

    if 'D0' in file:
        df = pd.DataFrame(
            columns=['Image_no', 'Trajectory', 'x_0', 'y_0', 'Radius'])
    elif 'D1' in file:
        df = pd.DataFrame(
            columns=['Image_no', 'Trajectory', 'x_1', 'y_1', 'Radius'])
    else:
        df = pd.DataFrame(
            columns=['Image_no', 'Trajectory', 'x', 'y', 'Radius'])
    return df


def save_df(file, cols, df):
    header = True
    df.to_csv('{}.csv'.format(file),
              mode='w',
              header=header,
              index_label=cols,
              index=False
              )


def concatenate_df(filename, sort_key='Trajectory'):
    # Check if filename-D0.csv and filename-D1.csv exist
    csv_files = list(
        map(lambda x: '{}{}.csv'.format(filename, x), ['-D0', '-D1']))

    if all(os.path.isfile(file) for file in csv_files):
        df0 = pd.read_csv(csv_files[0])
        df1 = pd.read_csv(csv_files[1])
        df = df0
        df0 = df0.set_index(['Image_no', 'Trajectory'])
        df1 = df1.set_index(['Image_no', 'Trajectory'])

        df = df.join(df1[['Y', 'Z']], on=[
                     'Image_no', 'Trajectory'], rsuffix='_D1')

        cols = df.columns.tolist()
        df = df[['Image_no', 'Trajectory', 'X', 'Y', 'Z', 'Z_D1', 'Radius']]
        if sort_key == 'Trajectory':
            df = df.sort_values(['Trajectory', 'Image_no'],
                                ascending=[True, True])
        elif sort_key == 'Image_no':
            df = df.sort_values(['Image_no', 'Trajectory'],
                                ascending=[True, True])
        df.to_csv('{}.csv'.format(filename),
                  mode='w',
                  header=True,
                  index_label=cols,
                  index=False
                  )


def open_overlaps_data(file):
    file = os.path.basename(file)
    return pd.read_csv(file + '_overlaps.csv')


def iterator(file, df, func, interpol_step, particle_no, *args):

    print('Analyzing file: {}'.format(file))
    averages = load_images_average(file)

    overlaps_data = open_overlaps_data(file)

    overlaps_data = overlaps_data.astype('int64')

    particles = overlaps_data['particle'].unique()

    for particle in particles:

        df = open_data_frame(file)
        cols = df.columns.tolist()
        trajectory_data = overlaps_data.loc[overlaps_data['particle'] ==
                                            particle, overlaps_data.columns != 'particle']
        trajectories = trajectory_data['trajectory'].tolist()
        for trajectory in trajectories:

            first_image = trajectory_data['start_frame'].loc[trajectory_data['trajectory'] ==
                                                             trajectory].iloc[0]
            last_image = trajectory_data['end_frame'].loc[trajectory_data['trajectory'] ==
                                                          trajectory].iloc[0]
            cx_real = trajectory_data['x_start'].loc[trajectory_data['trajectory'] ==
                                                     trajectory].iloc[0]
            cy_real = trajectory_data['y_start'].loc[trajectory_data['trajectory'] ==
                                                     trajectory].iloc[0]
            r = trajectory_data['r'].loc[trajectory_data['trajectory'] ==
                                         trajectory].iloc[0]

            print(first_image, last_image, cx_real, cy_real, r)
            images, ny, nx = load_images(file,
                                         first_image=first_image,
                                         last_image=last_image
                                         )

            start_time = time.time()
            for image_no, image in enumerate(images, first_image):
                # plt.imshow(image)
                # plt.show()
                T = Trajectory(cy_real, cx_real, r, image, averages,
                               ny, nx, image_no, trajectory, cols, r_real=r_real, step=interpol_step)

                cy_real, cx_real, df = func(T, df, *args)

                if (image_no - first_image) % 100 == 0:
                    print('{}\t{}\n'.format(image_no, time.time() - start_time))

        save_df(file + '_particle_{}'.format(particle), cols, df)


def second_pass(T, df):

    patch_average, bottom, left = average_image_of_patch(T)

    patch_average_2 = interpolate(patch_average, T.step)

    r_interpol = int(T.r * T.step)

    if T.image_no == 1:
        plt.imsave('{}{}.png'.format(filename.split(
            '/')[-1], T.cy_real), patch_average)

    cy, cx = center_of_darkest_circle(patch_average_2,
                                      r_interpol,
                                      interpol=True,
                                      cy=T.cy_real,
                                      cx=T.cx_real,
                                      bottom=bottom,
                                      left=left,
                                      step=T.step
                                      )

    cy_real = bottom + (-((T.step - 1) / 2) + cy) / T.step
    cx_real = left + (-((T.step - 1) / 2) + cx) / T.step

    row = [int(T.image_no), int(T.trajectory), cx_real, cy_real, T.r]
    dict_ = dict(zip(T.cols, row))
    df = df.append([dict_])

    return cy_real, cx_real, df


class Trajectory():
    def __init__(self, cy_real, cx_real, r, image, averages, ny, nx, image_no, trajectory, cols, step=1):
        self.cy_real = cy_real
        self.cx_real = cx_real
        self.r = r
        self.image = image
        self.averages = averages
        self.ny = ny
        self.nx = nx
        self.image_no = image_no
        self.trajectory = trajectory
        self.cols = cols
        self.step = step


if __name__ == "__main__":
    script, filename, interpol_step = argv

    main(filename,
         interpol_step=int(interpol_step))
