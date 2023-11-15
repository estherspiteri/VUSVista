import { IVus } from "../../models/view-vus.model.tsx/view-vus.model";

export interface IStoreAndVerifyVusFileRequest {
  vusFile: File;
}

export interface IStoreAndVerifyVusFileResponse {
  isSuccess: boolean;
}

export interface ILoadAllVusResponse {
  isSuccess: boolean;
  vusList?: IVus[] | null;
}
