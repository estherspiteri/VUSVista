import { IHomepageData } from "../../models/homepage.model";

export interface IHomepageResponse {
  isSuccess: boolean;
  data: IHomepageData;
}
